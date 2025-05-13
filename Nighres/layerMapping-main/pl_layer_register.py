'''
Mapping of quantitative T1 values from MP2RAGE data onto surface
=======================================

This is a pipeline processing BIDS formatted MP2RAGE data by performing the following steps:

01. Setup
02. MP2RAGE background cleaning
03. Resampling to 500 µm (hires pipeline only)
04. Skull stripping
05. Registration of whole brain to slab data (hires pipeline only)
06. Weighted image combination of whole brain and slab data (hires pipeline only)
07. FreeSurfer
08. Atlas-guided tissue classification using MGDM
09. Extract hemisphere
10. Crop volume to hemisphere
11. CRUISE cortical reconstruction
12. Extract layers across cortical sheet and map on surface
13. Process addtional data (if data is specified)
14. Process transform data (if data is specified)

Version 1.1. (23.02.2023) 
'''

############################################################################
# Contact information
# -------------------
# Dr. Falk Luesebrink
#
# falk dot luesebrink at ovgu dot de
# github.com/fluese
#

############################################################################
# Acknowledgement
# -------------------
# This work was supported by CRC 1436 "Neural Resources of Cognition" of
# the German Research Foundation (DFG) under project ID 425899996.

############################################################################
# THINGS TO DO
# -------------------
# 0. Largest Component for white segmentation?
# 1. Hard coded resolution for high resolution pipeline should be set automatically based on slab's resolution
# 2. Get rid of MATLAB:
#	   * MP2RAGE background cleaning [Needs to be re-written]
#	   * Combination of high resolution slab and low resolution whole brain data [Note: Function is included at the end. However, MATLAB function provides better (?) results. At least using the Python implementation CRUISE produced more error. Maybe due to using the FreeSurfer segmentation this is not an issue anymore]
# 3. Flag to write all data to disk or final results only
# 4. Add logging feature
#

############################################################################
# Installation instructions
# -------------------
# Things needed to be installed:
# 1. Nighres (tested with v1.4: https://nighres.readthedocs.io/en/latest/installation.html)
# 2. antspy (tested with v0.3.2: https://github.com/ANTsX/ANTsPy)
# 3. FreeSurfer (tested with v7.3.2: https://surfer.nmr.mgh.harvard.edu/fswiki/DownloadAndInstall)
# 4. MATLAB (tested with R2022a)
# 5. MP2RAGE-related-scripts (https://github.com/JosePMarques/MP2RAGE-related-scripts)
# 6. Tools for NIfTI and Analyze image format (https://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image)
# 7. Custom MATLAB scripts (weightedAverage and removeBackgroundnoise)
#

############################################################################
# Usage instructions
# -------------------
# This script requires whole brain MP2RAGE data organized according to BIDS.
# The script expects the first and second inversion, the T1 weighted data as
# well as the T1 map. You can make use of a high resolution MP2RAGE slab
# additionally. The slab will be merged into the whole brain volume
# automatically and used for further processing. To make use of a high 
# resolution slab, the slab needs to be acquired as the second run in a
# session (or at least named as if it were acquired in the same session).
#
# In the subsection "set parameters" of the setup section below, you need to
# specify the folder to the BIDS directory and the label of the subject you
# want to process. This can be done in multiple ways. First, specify the 
# path and file in the script directly. Secondly, when the script is called
# the next two inputs are for the subject ID and path to the BIDS directory,
# respectively. Lastly, if all is left empty, the user will be asked to
# specify a subject ID and the full path to the BIDS directory at the 
# beginning of the script.
#
# Example use (Case 1):
# python3 layerMapping.py aaa /path/to/BIDS/data/
#
# Example use (Case 2):
# python3 layerMapping.py aaa
#
# Example use (Case 3):
# python3 layerMapping.py
#
# In case 1 the subject with the ID 'aaa' will be processed. The according
# data is expected in '/path/to/BIDS/data/'. Output will be generated in
# '/path/to/BIDS/data/derivatives/sub-aaa/'. In case 2, either the path
# to the data needs to be defined in the script itself in the variable
# 'BIDS_path' or if left empty the user is asked to specify the path
# to the data. In case 3, the same for the path is true for the subject
# ID. It can either be specified in the script itself in the 'sub'
# variable or if left empty, the user will be asked for an ID during
# processing.
#
# Furthermore, you can specify a single volume which will be mapped on the
# surface additionally. This is meant for mapping another contrast  on the
# surface, e.g. QSM data or fMRI results. The resulting transformation can
# be applied either to another volume or to all compressed NIfTI files of
# an entire directory. Currently, this can be done in the script only.
#
# To enable logging of the output, one can send the screen output into a
# file. This can be done, e.g. by using following command:
# python3 -u layerMapping.py aaa |& tee -a /tmp/luesebrink/derivatives/sub-aaa/layerMapping.log
#
# A few flags have been set up which allow to reprocess data at various
# stages of the pipeline (e.g. entirely, from the segmentation onwards,
# mapping of additional data, or per hemisphere). In case you want to 
# map multiple files onto the surface additionally, this is especially 
# useful as you don't have to process all other data again. Simply
# specify a file under 'map_data' that shall be mapped onto the surface
# and change the 'reprocess_map_data' flag to True.
#

############################################################################
# 1. Setup
# -------------------
# First, we import nighres and all other modules. Then, we will define all
# required parameters.

############################################################################
# 1.1. Import python libraries
# -------------------
#import nighres
import os
import sys
import ants
import nibabel as nb
import glob
import subprocess
from nilearn.image import mean_img
from nilearn.image import crop_img
from nibabel.processing import conform
from time import localtime, strftime

############################################################################
# 1.2. Set parameters
# -------------------
# Define subject to be processed according to BIDS nomenclature here. If
# left empty, the user will be asked to specify the subject's ID during
# processing. Otherwise, the first input after the script is meant to
# specify the subject ID.
sub = 'unk'

# Define BIDS path here. If left empty, the user will be asked to specify
# the path during processing. Otherwise, the second input after the script
# is meant to specify the BIDS_path. In a future release this will be 
# changed for easier handling.
BIDS_path = '/Volumes/IKND/AG_Kuehn/Peng/LayerPRF/LayerPreparation/SurfaceRegister/unk742/'

# Process with an additional high resolution MP2RAGE slab. If 'True' the 
# first run must be the lower resolution whole brain MP2RAGE volume and
# the second run must be the higher resolution MP2RAGE slab volume.
hires = False

# Map specific volume onto the surface. This could be the BOLD of a task
# fMRI time series (preferably the mean across the time series) or the
# magnitude data of a QSM volume. This volume will be registered to the
# T1map of the MP2RAGE volume.
# Requires an absolute path to a NIfTI file. Results will be written to
# <BIDS_path>/derivatives/sub-<label>/. If the path points to a 
# non-existing file, the according option will be omitted.
map_data = BIDS_path + 'derivatives/sub-unk/adata_mean.nii.gz'

# Transform data in the same space as 'map_data' which is then mapped onto
# the surface. Could for example be the statistical maps of SPM from fMRI or the
# Chi map of QSM data.
# Requires an absolute path to a NIfTI file or a directory containing compressed
# NIfTI files. The transformation is then applied to all files within the
# directory. Results will be written to <BIDS_path>/derivatives/sub-<label>/.
# If the path points to a non-existing file or directory, the according
# option will be omitted.
transform_data = BIDS_path + 'derivatives/sub-unk/Fourier_Mean/'

# Choose interpolation method for mapping of additional data. Choice of
# interpolator can be 'linear', 'nearstNeighbour, 'bSpline', 'genericLabel'.
# See https://antspy.readthedocs.io/en/latest/_modules/ants/registration/apply_transforms.html 
# for full list of supported interpolators.
interpolation_method = 'nearestNeighbor'

# Flag to overwrite all existing data
reprocess = False

# Flag to start reprocessing from FreeSurfer onwards
reprocess_FreeSurfer = False

# Flag to start reprocessing from segmentation onwards
reprocess_segmentation = False

# Flag to reprocess left or right hemisphere only
reprocess_leftHemisphere = False
reprocess_rightHemisphere = False

# Flag to reprocess additional mapping (and transformation) data. This flag
# is especially useful if you want to map additional information onto the
# surface without re-running the entire pipeline.
# This could be the case if you have mapped functional data on the surface
# already already and now you also want to map QSM data onto the surface. 
# In this scenario you would want to set this flag to true, change the path
# of 'map_data' and 'transform_data' to your QSM data and the resulting
# Chi map, respectively. All other reprocess flags should be set to 'False'
# to avoid unnessecary reprocessing.
# The basename of the output will be based on the file name of the input
# data. According to BIDS the subject label will be added as prefix to all
# output data.
reprocess_map_data = False

# Define path from where data is to be copied into "BIDS_path". Could either
# be used to create a backup of the data before processing or transfer to a
# compute server.
# If the path does not exist, this option will be omitted and data from
# "BIDS_path" will be used instead.
# Unfortunately, does not work for network drives. -> Workaround?
copy_data_from = ''

# Define path to atlas. Here, we use custom priors which seem to work well
# with 7T MP2RAGE data collected in Magdeburg, Germany. Otherwise, please 
# choose the recent priors that come along with nighres (currently 3.0.9)
# instead of using the default atlas (currently 3.0.3). You can find that
# text file in the python package of nighres in the atlas folder.
atlas = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'atlas', 'brain-atlas-quant-3.0.9_customPriors.txt')

############################################################################
# 1.3. Define file names (and optionally create backup or copy data to
# compute server)
# -------------------
if not sub:
	if len(sys.argv) > 1:
		sub = sys.argv[1]
	else:
		sub = input('Please enter subject ID: ')
	
if not BIDS_path:
	if len(sys.argv) > 2:
		BIDS_path = sys.argv[2]
	else:
		BIDS_path = input('Please enter full path to BIDS directory: ')

# Set paths and create directories
in_dir = BIDS_path + 'sub-' + sub + '/anat/'
out_dir = BIDS_path + 'derivatives/sub-' + sub + '/'
subprocess.run(["mkdir", "-p", in_dir])
subprocess.run(["mkdir", "-p", out_dir])

print('*****************************************************')
print('* Processing pipeline started at: ' + strftime("%Y-%m-%d %H:%M:%S", localtime()))
print('*****************************************************')
print('Working directory: ' + BIDS_path)
print('Subject to be processed: ' + sub)
open(os.path.join(out_dir, 'is_running'), mode='a').close()

# Define file names
if hires == True:
	INV1 = os.path.join(in_dir, 'sub-' + sub + '_run-01_inv-1_MP2RAGE.nii.gz')
	INV2 = os.path.join(in_dir, 'sub-' + sub + '_run-01_inv-2_MP2RAGE.nii.gz')
	T1map = os.path.join(in_dir, 'sub-' + sub + '_run-01_T1map.nii.gz')
	UNI = os.path.join(in_dir, 'sub-' + sub + '_run-01_UNIT1.nii.gz')

	INV1_slab = os.path.join(in_dir, 'sub-' + sub + '_run-02_inv-1_MP2RAGE.nii.gz')
	INV2_slab = os.path.join(in_dir, 'sub-' + sub + '_run-02_inv-2_MP2RAGE.nii.gz')
	T1map_slab = os.path.join(in_dir, 'sub-' + sub + '_run-02_T1map.nii.gz')
	UNI_slab = os.path.join(in_dir, 'sub-' + sub + '_run-02_UNIT1.nii.gz')
else:
	INV1 = os.path.join(in_dir, 'sub-' + sub + '_run-01_inv-1_MP2RAGE.nii.gz')
	INV2 = os.path.join(in_dir, 'sub-' + sub + '_run-01_inv-2_MP2RAGE.nii.gz')
	T1map = os.path.join(in_dir, 'sub-' + sub + '_run-01_T1map.nii.gz')
	UNI = os.path.join(in_dir, 'sub-' + sub + '_run-01_UNIT1.nii.gz')

# Copy data if path exists
print('')
print('*****************************************************')
print('* Data transfer to working directory.')
print('* Started at: ' + strftime("%Y-%m-%d %H:%M:%S", localtime()))
print('*****************************************************')
if os.path.isdir(copy_data_from):
	if os.path.isfile(UNI)  and reprocess != True:
		print('Files exists already. Skipping data transfer.')
	else:
		subprocess.run(["scp", "-r", copy_data_from, " ", BIDS_path])
else:
	print('WARNING: No path found to copy data from found. Continuing without copying data!')

# Check if data exists.
print('')
print('*****************************************************')
print('* Checking if files exists.')
print('* Started at: ' + strftime("%Y-%m-%d %H:%M:%S", localtime()))
print('*****************************************************')
if os.path.isfile(INV1) and os.path.isfile(INV2) and os.path.isfile(T1map) and os.path.isfile(UNI):
	print('Files exists. Good!')
else:
	print('No files found. Exiting!')
	print('')
	print('INFO: Be sure that your input data are structured according to BIDS.')
	print('      It should be organized like this:')
	print('')
	print('		 ./BIDS_path/')
	print(' 	 └── /sub-' + sub + '/')
	print('		 	  └── anat')
	print('			      ├── sub-' + sub + '_run-01_inv-1_MP2RAGE.nii.gz')
	print('			      ├── sub-' + sub + '_run-01_inv-2_MP2RAGE.nii.gz')
	print('			      ├── sub-' + sub + '_run-01_T1map.nii.gz')
	print('			      ├── sub-' + sub + '_run-01_UNIT1.nii.gz')
	print('			      ├── sub-' + sub + '_run-02_inv-1_MP2RAGE.nii.gz')
	print('			      ├── sub-' + sub + '_run-02_inv-2_MP2RAGE.nii.gz')
	print('			      ├── sub-' + sub + '_run-02_T1map.nii.gz')
	print('			      └── sub-' + sub + '_run-02_UNIT1.nii.gz')
	print('')
	exit()

if hires == True:
	if os.path.isfile(INV1_slab) and os.path.isfile(INV2_slab) and os.path.isfile(T1map_slab) and os.path.isfile(UNI_slab):
		print('High resolution files exists. Good!')
	else:
		print('WARNING: No high resolution files found. Continuing without hires option!')
		hires = False

# Check if file for mapping of additional data exists.
print('')
print('*****************************************************')
print('* Checking for additional data to be mapped.')
print('* Started at: ' + strftime("%Y-%m-%d %H:%M:%S", localtime()))
print('*****************************************************')
if os.path.isfile(map_data):
	map_file_onto_surface = True
	print('File found for mapping of additional data onto surface. Good!')

	# Change name string for output of map data. This will first remove the path
	# of given the file for an appropriate naming of the output data.
	file_name = os.path.basename(map_data)
	index_of_dot = file_name.index('.')
	map_data_output = file_name[:index_of_dot]
	# This will remove the sub-<label> string if the data is already BIDS
	# formatted. It'll be added during creation of the output data later
	# again.
	if 'sub-' + sub in map_data_output:
		map_data_output = map_data_output[8:]
else:
	map_file_onto_surface = False
	print('WARNING: Could not find file for surface mapping. Omitting flag!')

# Check if file for applying the transformation of the additional data exists.
if os.path.isfile(transform_data) == False and os.path.isdir(transform_data) == False:
	map_transform_file_onto_surface = False
	print('WARNING: Could not find file or path for applying transform of additional data. Omitting flag!')
elif os.path.isdir(transform_data) or os.path.isfile(transform_data):
	map_transform_file_onto_surface = True
	print('File or path found for applying transform of additional data. Good!')

	# Change name string for output of transformed data. This will remove the 
	# path of the input given for an appropriate naming of the output data.
	# In case a path is given as input a list with all compressed NIfTI files
	# will be created. Then the path from that list will be stripped to create
	# a list for appropriate naming of the output data.
	if os.path.isdir(transform_data):
		transform_data = glob.glob(transform_data + '*.nii.gz')
		transform_data_output = []
		if isinstance(transform_data, list):
			for tmp in transform_data:
				file_name = os.path.basename(tmp)
				index_of_dot = file_name.index('.')
				tmp = file_name[:index_of_dot]
				# This will remove the sub-<label> string if the data is already BIDS
				# formatted. It'll be added during creation of the output data later
				# again.
				if 'sub-' + sub in tmp:
					tmp = tmp[8:]
				transform_data_output.append(tmp)
	else:
		file_name = os.path.basename(transform_data)
		index_of_dot = file_name.index('.')
		transform_data_output = file_name[:index_of_dot]
		# This will remove the sub-<label> string if the data is already BIDS
		# formatted. It'll be added during creation of the output data later
		# again.	
		if 'sub-' + sub in transform_data_output:
			transform_data_output = transform_data_output[8:]

# For naming and checking of processing files.
if hires == True:
	merged = '_merged_run-01+02'
	resampled = '_resampled'
else:
	merged = '_run-01'
	resampled = ''

############################################################################
# 2. MP2RAGE background cleaning
# -------------------
# This script creates MP2RAGE T1w images without the strong background noise in
# air regions as implemented by Marques [Taken from his Github repository] using
# the method of O'Brien.
#
# O'Brien, et al, 2014.
# Robust T1-Weighted Structural Brain Imaging and Morphometry at 7T Using MP2RAGE
# PLOS ONE 9, e99676. doi:10.1371/journal.pone.0099676
# http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0099676
# -------------------
print('')
print('*****************************************************')
print('* MP2RAGE background cleaning')
print('* Started at: ' + strftime("%Y-%m-%d %H:%M:%S", localtime()))
print('*****************************************************')
# Check for file and start MATLAB to conduct background cleaning.
if os.path.isfile(os.path.join(out_dir, 'sub-' + sub + '_run-01_T1w.nii.gz')) and reprocess != True:
	print('File exists already. Skipping process.')
else:
	reg = str(10)
	output = os.path.join(out_dir, 'sub-' + sub + '_run-01_T1w.nii.gz')
	cmd = "removeBackgroundnoise(\'" + UNI + "\', \'" + INV1 + "', \'" + INV2 + "\', \'" + output + "\'," + reg + "); exit;"
	subprocess.run(["matlab", "-nosplash", "-nodisplay", "-r", cmd])

# Check for high resolution slab data and start MATLAB to conduct background cleaning.
if hires == True:
	if os.path.isfile(os.path.join(out_dir, 'sub-' + sub + '_run-02_T1w.nii.gz')) and reprocess != True:
		print('File exists already. Skipping process.')
	else:
		output = os.path.join(out_dir, 'sub-' + sub + '_run-02_T1w.nii.gz')
		cmd = "removeBackgroundnoise(\'" + UNI_slab + "\', \'" + INV1_slab + "', \'" + INV2_slab + "\', \'" + output + "\'," + reg + "); exit;"
		subprocess.run(["matlab", "-nosplash", "-nodisplay", "-r", cmd])
		
# Update file names
T1w = os.path.join(out_dir, 'sub-' + sub + '_run-01_T1w.nii.gz')
T1w_slab = os.path.join(out_dir, 'sub-' + sub + '_run-02_T1w.nii.gz')

############################################################################
# 3. Resampling to higher resolution data
# -------------------
# Resample image to an isotropic resolution of 500 µm.
#
# Should be changed to the resolution of the high resolution data at some
# point and not expect 500 µm data.
if hires == True:
	print('')
	print('*****************************************************')
	print('* Resampling MP2RAGE to an isotropic resolution of 500 µm')
	print('* Started at: ' + strftime("%Y-%m-%d %H:%M:%S", localtime()))
	print('*****************************************************')

	# Check for file and resample T1w and T1map using 4th order b spline interpolation.
	if os.path.isfile(os.path.join(out_dir, 'sub-' + sub + '_run-01_T1map_resampled.nii.gz')) and reprocess != True:
		print('File exists already. Skipping process.')
	else:
		input_img = nb.load(T1w)
		resampled_img = conform(input_img, out_shape=(336,448,448), voxel_size=(0.5,0.5,0.5), order=4)
		nb.save(resampled_img, os.path.join(out_dir, 'sub-' + sub + '_run-01_T1w_resampled.nii.gz'))

		input_img = nb.load(T1map)
		resampled_img = conform(input_img, out_shape=(336,448,448), voxel_size=(0.5,0.5,0.5), order=4)
		nb.save(resampled_img, os.path.join(out_dir, 'sub-' + sub + '_run-01_T1map_resampled.nii.gz'))
		print('Done.')

	# Update file names
	T1map = os.path.join(out_dir, 'sub-' + sub + '_run-01_T1map_resampled.nii.gz')
	T1w = os.path.join(out_dir, 'sub-' + sub + '_run-01_T1w_resampled.nii.gz')
else:
	# Copy and update file name
	subprocess.run(["cp", T1map, os.path.join(out_dir, "sub-" + sub + "_run-01_T1map.nii.gz")])
	T1map = os.path.join(out_dir, 'sub-' + sub + '_run-01_T1map.nii.gz')

############################################################################
# 4. Skull stripping
# ----------------
# Here, we perform skull stripping of T1w data using FreeSurfer's
# mri_synthstrip routine. Afterwards, we apply the mask to the T1map
# data for consistency.
print('')
print('*****************************************************')
print('* Skull stripping of whole brain data')
print('* Started at: ' + strftime("%Y-%m-%d %H:%M:%S", localtime()))
print('*****************************************************')
T1w = os.path.join(out_dir, 'sub-' + sub + '_run-01_T1w' + resampled + '.nii.gz')
T1w_masked = os.path.join(out_dir, 'sub-' + sub + '_run-01_T1w' + resampled + '_masked.nii.gz')
mask = os.path.join(out_dir, 'sub-' + sub + '_run-01_T1w' + resampled + '_brainmask.nii.gz')

if os.path.isfile(T1w_masked) and reprocess != True:
	print('File exists already. Skipping process.')
else:
	subprocess.run(["mri_synthstrip", "-i", T1w, "-o", T1w_masked, "-m", mask, "-b", "-1"])

if hires == True:
	print('')
	print('*****************************************************')
	print('* Skull stripping of high resolution slab data')
	print('* Started at: ' + strftime("%Y-%m-%d %H:%M:%S", localtime()))
	print('*****************************************************')
	T1w_slab = os.path.join(out_dir, 'sub-' + sub + '_run-02_T1w.nii.gz')
	T1w_slab_masked = os.path.join(out_dir, 'sub-' + sub + '_run-02_T1w_masked.nii.gz')
	mask_slab = os.path.join(out_dir, 'sub-' + sub + '_run-02_T1w_brainmask.nii.gz')
	
	if os.path.isfile(T1w_slab_masked) and reprocess != True:
		print('File exists already. Skipping process.')
	else:
		subprocess.run(["mri_synthstrip", "-i", T1w_slab, "-o", T1w_slab_masked, "-m", mask_slab, "-b", "-1"])

print('')
print('*****************************************************')
print('* Masking T1map of whole brain data')
print('* Started at: ' + strftime("%Y-%m-%d %H:%M:%S", localtime()))
print('*****************************************************')
# Update file names	
T1w = os.path.join(out_dir, 'sub-' + sub + '_run-01_T1w' + resampled + '.nii.gz')
T1w_masked = os.path.join(out_dir, 'sub-' + sub + '_run-01_T1w' + resampled + '_masked.nii.gz')
T1map = os.path.join(out_dir, 'sub-' + sub + '_run-01_T1map' + resampled + '.nii.gz')
T1map_masked = os.path.join(out_dir, 'sub-' + sub + '_run-01_T1map' + resampled + '_masked.nii.gz')
mask = os.path.join(out_dir, 'sub-' + sub + '_run-01_T1w' + resampled + '_brainmask.nii.gz')

# Check if volume exists
if os.path.isfile(T1map_masked) and reprocess != True:
	print('File exists already. Skipping process.')
else:
# Mask image otherwise
	maskedImage = ants.mask_image(ants.image_read(T1map), ants.image_read(mask))
	ants.image_write(maskedImage, os.path.join(out_dir, 'sub-' + sub + '_run-01_T1map' + resampled + '_masked.nii.gz'))
	print('Done.')

if hires == True:
	print('')
	print('*****************************************************')
	print('* Masking T1map of high resolution slab data')
	print('* Started at: ' + strftime("%Y-%m-%d %H:%M:%S", localtime()))
	print('*****************************************************')
	# Update file names	
	T1w_slab = os.path.join(out_dir, 'sub-' + sub + '_run-02_T1w.nii.gz')
	T1w_slab_masked = os.path.join(out_dir, 'sub-' + sub + '_run-02_T1w_masked.nii.gz')
	T1map_slab = os.path.join(in_dir, 'sub-' + sub + '_run-02_T1map.nii.gz')
	T1map_slab_masked = os.path.join(out_dir, 'sub-' + sub + '_run-02_T1map_masked.nii.gz')
	mask_slab = os.path.join(out_dir, 'sub-' + sub + '_run-02_T1w_brainmask.nii.gz')

	# Check if volume exists
	if os.path.isfile(T1map_slab_masked) and reprocess != True:
		print('File exists already. Skipping process.')
	else:
	# Mask image otherwise
		subprocess.run(["cp", os.path.join(in_dir, "sub-" + sub + "_run-02_T1map.nii.gz"), os.path.join(out_dir, "sub-" + sub + "_run-02_T1map.nii.gz")])
		maskedImage = ants.mask_image(ants.image_read(T1map_slab), ants.image_read(mask_slab))
		ants.image_write(maskedImage, os.path.join(out_dir, 'sub-' + sub + '_run-02_T1map_masked.nii.gz'))
		print('Done.')

############################################################################
# 13. Process additional data.
# -------------------
if reprocess_map_data:
  reprocess = True

############################################################################
# 13.1. Register additional data to (upsampled) T1map
# -------------------
if map_file_onto_surface:
	if hires == True:
		print('')
		print('*****************************************************')
		print('* Register additional data to upsampled T1map')
		print('* Started at: ' + strftime("%Y-%m-%d %H:%M:%S", localtime()))
		print('*****************************************************')
		reg1 = 'sub-' + sub + '_' + map_data_output + '_registered_to_sub-' + sub + '_run-01_T1map_resampled'
	else:
		print('')
		print('*****************************************************')
		print('* Register additional data to T1map')
		print('* Started at: ' + strftime("%Y-%m-%d %H:%M:%S", localtime()))
		print('*****************************************************')
		reg1 = 'sub-' + sub + '_' + map_data_output + '_registered_to_sub-' + sub + '_run-01_T1map'

	if os.path.isfile(os.path.join(out_dir, reg1 + '.nii.gz')) and reprocess != True:
		print('File exists already. Skipping process.')
		map_data = os.path.join(out_dir, reg1 + '.nii.gz')
	else:
		# Here, we mask the data to be registered to improve the registration
		map_data_mask = ants.image_read(map_data)
		map_data_mask = ants.get_mask(map_data_mask)
		map_data_masked = ants.mask_image(ants.image_read(map_data),map_data_mask)
		ants.image_write(map_data_masked, out_dir + 'sub-' + sub + '_' + map_data_output + '_masked.nii.gz')

		# Register additional data non-linearly to T1map using mutual information as similarity metric
		registeredImage = ants.registration(
			fixed = map_data_masked,
			moving = ants.image_read(T1map_masked),
			type_of_transform = 'SyNRA',
			aff_random_sampling_rate = 1,
			grad_step = 0.1,
			reg_iterations = (200, 100, 50 , 30),
			verbose = True,
			# outprefix = out_dir + reg1 + '_',
			)

		# Apply transformation
		warpedImage = ants.apply_transforms(
			fixed = ants.image_read(T1map),
			moving = ants.image_read(map_data),
			transformlist = registeredImage['invtransforms'],
			interpolator = interpolation_method,
			verbose = True,
			)

		# Write file to disk
		ants.image_write(warpedImage, out_dir + reg1 + '.nii.gz') 

		# Update file name
		map_data = os.path.join(out_dir, reg1 + '.nii.gz')

############################################################################
# 14. Process transformed data.
# -------------------
if reprocess_map_data:
	reprocess = True

############################################################################
# 14.1 Apply transformation to data to be mapped on the surface
# -------------------
# Define file name for output of transformation.
if map_file_onto_surface and map_transform_file_onto_surface:
	if isinstance(transform_data_output, list):
		reg2 = []
		for name in transform_data_output:
			if hires == True:
				reg2.append('sub-' + sub + '_' + name + '_registered_to_' + sub + '_run-01_T1map_resampled.nii.gz')
			else:
				reg2.append('sub-' + sub + '_' + name + '_registered_to_' + sub + '_run-01_T1map.nii.gz')

	else:
		if hires == True:
			reg2 = 'sub-' + sub + '_' + transform_data_output + '_registered_to_' + sub + '_run-01_T1map_resampled.nii.gz'
		else:
			reg2 = 'sub-' + sub + '_' + transform_data_output + '_registered_to_' + sub + '_run-01_T1map.nii.gz'

	print('')
	print('*****************************************************')
	print('* Apply transformation of registration to additional data')
	print('* Started at: ' + strftime("%Y-%m-%d %H:%M:%S", localtime()))
	print('*****************************************************')
	# Check if a list exists and if so, if the last file of list exists already.
	if isinstance(transform_data, list):
		if os.path.isfile(os.path.join(out_dir, reg2[-1])) and reprocess != True:
			print('File exists already. Skipping process.')
		else:
			# Iterate over length of list
			length = len(transform_data)
			for i in range(length):
				# Apply transformation
				warpedImage = ants.apply_transforms(
					fixed = ants.image_read(T1map),
					moving = ants.image_read(transform_data[i]),
					transformlist = registeredImage['invtransforms'],
					interpolator = interpolation_method,
					#imagetype = 3,
					verbose = True,
					)

				# Write file to disk
				ants.image_write(warpedImage, out_dir + reg2[i])

	else:
	# Check if output exists already.
		if os.path.isfile(os.path.join(out_dir, reg2)) and reprocess != True:
			print('File exists already. Skipping process.')
		else:
			# Apply transformation
			warpedImage = ants.apply_transforms(
				fixed = ants.image_read(T1map),
				moving = ants.image_read(transform_data),
			    transformlist = registeredImage['invtransforms'],
				interpolator = interpolation_method,
				verbose = True,
				)

			# Write file to disk
			ants.image_write(warpedImage, out_dir + reg2) 
		
	# Update file name
	if isinstance(reg2, list):
		transform_data = []
		for tmp in reg2:
			transform_data.append(os.path.join(out_dir, tmp))		
	else:
		transform_data = os.path.join(out_dir, reg2)
