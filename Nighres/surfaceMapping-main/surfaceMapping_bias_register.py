'''
Mapping of quantitative T1 values from MP2RAGE data onto surface
=======================================

This is a pipeline processing MP2RAGE data by performing the following steps:

01. Setup
02. MP2RAGE background cleaning
03. Resampling to 500 µm (only for hires option)
04. Imhomogeneity correction and skull stripping
05. Registration of whole brain to slab data (optionally)
06. Regisration of additional to structural data (optionally)
07. Weighted image combination of whole brain and slab data (optionally)
08. Atlas-guided tissue classification using MGDM
09. Region extraction (left hemisphere) 
10. Crop volume (left hemisphere)
11. CRUISE cortical reconstruction (left hemisphere)
12. Extract layers across cortical sheet and map on surface (left hemisphere)
13. Region extraction (right hemisphere)
14. Crop volume (right hemisphere)
15. CRUISE cortical reconstruction (right hemisphere)
16. Extract layers across cortical sheet and map on surface (right hemisphere)

Version 0.96 (03.09.2022)
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
# THINGS TO DO
# -------------------
# 1. Specify files, paths, and flag from outside this script
# 2. Get rid of MATLAB:
# 	Switch from SPM's bias correction method to N4 of ANTs for easier use
# 	Re-write the weightedAverage script for Python
#	MP2RAGE background cleaning
# 3. Flag to write all data to disk or final results only
# 4. Add logging feature
#

############################################################################
# Installation instructions
# -------------------
# Things needed to be installed:
# 1. Nighres (https://nighres.readthedocs.io/en/latest/installation.html)
# 2. antspy (https://github.com/ANTsX/ANTsPy)
# 3. MATLAB
# 4. SPM12 (https://www.fil.ion.ucl.ac.uk/spm/software/download/)
# 5. MP2RAGE-related-scripts (https://github.com/JosePMarques/MP2RAGE-related-scripts)
# 6. Custom MATLAB scripts (weightedAverage, removeBackgroundnoise and biasCorrection)
#
# You need to change the path to the tissue probability model for the bias
# field correction method. This needs to be done in
# 
# ./biasCorrection/preproc_sensemap.m on line 19
#

############################################################################
# Usage instructions
# -------------------
# This script requires whole brain MP2RAGE data organized according to BIDS.
# You can make use of a high resolution MP2RAGE slab additionally. The slab
# will be merged into the whole brain volume automatically and used for
# further processing. 
#
# In the subsection "set parameters" of the setup section below, you need to
# specify the folder to the BIDS directory and the label of the subject you
# want to process.
#
# Furthermore, you can specify a single volume which will be mapped on the
# surface additionally. This is meant for mapping another contrast  on the
# surface, e.g. QSM data or fMRI results. The resulting transformation can
# be applied either to another volume or to all compressed NIfTI files of
# an entire directory.
#
# To enable logging of the output, one can send the screen output into a
# file. This can be done, e.g. by using following command:
# python3 -u surfaceMapping.py |& tee -a /tmp/luesebrink/surfaceMapping.log
#
# A few flags have been set up which allow to reprocess data at various
# stages of the pipeline (e.g. entirely, from the segmentation onwards,
# mapping of additional data, or per hemisphere).
#

############################################################################
# 1. Setup
# -------------------
# First, we import nighres, the os and nibabel module. Set the in- and output
# directory and define the file names. In the final release this will be
# changed accordingly.

############################################################################
# 1.1. Import python libraries
# -------------------
#import nighres
import os
import sys
import ants
import nibabel as nb
import glob
from nilearn.image import mean_img
from nilearn.image import crop_img
from nibabel.processing import conform
from time import localtime, strftime

############################################################################
# 1.2. Set parameters
# -------------------
# Define BIDS path
BIDS_path = '/Users/pliu/Documents/LayerPRF/SurfaceMapping/hby152/'

# Define subject following BIDS
sub = 'hby'

# Process with an additional high resolution MP2RAGE slab. If 'True' the 
# first run must be the lower resolution whole brain MP2RAGE volume and
# the second run must be the higher resolution MP2RAGE slab volume.
hires = True

# Map specific volume onto the surface. This could be the BOLD of a task
# fMRI time series (preferably the mean across the time series) or the
# magnitude data of a QSM volume. This volume will be registered to the
# T1map of the MP2RAGE volume.
# Requries an absolute path to a NIfTI file. Results will be written to
# <BIDS_path>/derivatives/sub-<label>/. If the path points to a 
# non-existing file, the according option will be omitted.
#map_data = BIDS_path + 'derivatives/sub-qxo/func/sub-qxo_task-pRF_bold_D3_mean.nii.gz'
map_data = BIDS_path + 'derivatives/sub-hby/func/sub-hby_task-pRF_bold_mean.nii.gz'
#map_data = BIDS_path + 'derivatives/sub-qxo/func/sub-qxo_task-Loc_bold_mean.nii.gz'
#map_data = ''

# Transform data in the same space as 'map_data' which is then mapped onto
# the surface. Could for example be the statistical maps of SPM from fMRI or QSM
# data.
# Requries an absolute path to a NIfTI file or a directory containing compressed
# NIfTI files. The transformation is then applied to all files within the
# directory. Results will be written to <BIDS_path>/derivatives/sub-<label>/.
# If the path points to a non-existing file or directory, the according
# option will be omitted.
#transform_data = BIDS_path + 'derivatives/sub-qxo/pRF/pRF_D3/'
#transform_data = BIDS_path + 'derivatives/sub-qxo/pRF/Phase_D2+D3/'
#transform_data = BIDS_path + 'derivatives/sub-qxo/Loc/'
transform_data = ''

# Choose interpolation method for mapping of additional data. Choice of interpolator can be 'linear', 'nearstNeighbor, 'bSpline', 'genericLabel'. See https://antspy.readthedocs.io/en/latest/_modules/ants/registration/apply_transforms.html for full list of supported interpolators.
interpolation_method = 'nearestNeighbor'

# Flag to overwrite all existing data
reprocess = False

# Flag to start reprocessing from segmentation onwards
reprocess_segmentation = False

# Flag to reprocess left or right hemisphere only
reprocess_leftHemisphere = False
reprocess_rightHemisphere = False

# Flag to reprocess additional mapping (and transformation) data. This flag
# is especially useful if you want to map more information onto the surface
# without re-running the entire pipeline. The basename of the output will be
# based on the file name of the input data. According to BIDS the subject
# label will be added as prefix to all output data.
reprocess_map_data = False

# Define path from where data is to be copied into "BIDS_path". Could either
# be to create a backup of the data before processing or transfer to a
# compute server.
# If the path does not exist, this option will be omitted and the data from
# "BIDS_path" will be used instead. Does not work for network drives.
copy_data_from = 'gerd:/media/luesebrink/bmmr_data/data/sensemap/sample_data/'

# Define path to atlas. Here, we use custom priors which seem to work well
# with 7T MP2RAGE data collected in Magdeburg, Germany. Otherwise, please 
# choose the recent priors that come along with nighres (currently 3.0.9)
# instead of using the default atlas (currently 3.0.3). You can find that
# text file in the python package of nighres in the atlas folder.
atlas = '/surfacemapping/atlas/brain-atlas-quant-3.0.9_customPriors.txt'

############################################################################
# 1.2.1. Display stuff:
# -------------------
print('*****************************************************')
print('* Processing pipeline started at: ' + strftime("%Y-%m-%d %H:%M:%S", localtime()))
print('*****************************************************')
print('Working directory: ' + BIDS_path)
print('Subject to be processed: ' + sub)

############################################################################
# 1.3. Define file names (and optionally create backup or copy data to
# compute server)
# -------------------
# Set paths
in_dir = BIDS_path + 'sub-' + sub + '/anat/'
out_dir = BIDS_path + 'derivatives/sub-' + sub + '/'
os.system('mkdir -p ' + in_dir)
os.system('mkdir -p ' + out_dir)

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
if os.path.isdir(copy_data_from):
	print('')
	print('*****************************************************')
	print('* Data transfer to working directory.')
	print('* Started at: ' + strftime("%Y-%m-%d %H:%M:%S", localtime()))
	print('*****************************************************')
	if os.path.isfile(UNI)  and reprocess != True:
		print('Files exists already. Skipping data transfer.')
	else:
		os.system('scp -r ' + copy_data_from + ' ' + BIDS_path)
else:
	print('')
	print('*****************************************************')
	print('* Data transfer to working directory.')
	print('* Started at: ' + strftime("%Y-%m-%d %H:%M:%S", localtime()))
	print('*****************************************************')
	print('No path found to copy data from found. Continuing without copying data!')

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
	exit()

if hires == True:
	if os.path.isfile(INV1_slab) and os.path.isfile(INV2_slab) and os.path.isfile(T1map_slab) and os.path.isfile(UNI_slab):
		print('High resolution files exists. Good!')
	else:
		print('No high resolution files found. Continuing without hires option!')
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
	print('Could not find file for surface mapping. Omitting flag!')

# Check if file for applying the transformation of the additional data exists.
if os.path.isfile(transform_data) == False and os.path.isdir(transform_data) == False:
	map_transform_file_onto_surface = False
	print('Could not find file or path for applying transform of additional data. Omitting flag!')
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

# Set flag if high resolution slab is used.
if hires == True:
	filename = '_merged_run-01+02'
else:
	filename = '_run-01'

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
	os.system('matlab -nosplash -nodisplay -r \"removeBackgroundnoise(\'' + UNI + '\', \'' + INV1 + '\', \'' + INV2 + '\', \'' + output + '\', ' + reg + '); exit;\"')

# Check for high resolution slab data and start MATLAB to conduct background cleaning.
if hires == True:
	if os.path.isfile(os.path.join(out_dir, 'sub-' + sub + '_run-02_T1w.nii.gz')) and reprocess != True:
		print('File exists already. Skipping process.')
	else:
		output = os.path.join(out_dir, 'sub-' + sub + '_run-02_T1w.nii.gz')
		os.system('matlab -nosplash -nodisplay -r \"removeBackgroundnoise(\'' + UNI_slab + '\', \'' + INV1_slab + '\', \'' + INV2_slab + '\', \'' + output + '\', ' + reg + '); exit;\"')

# Update file names
T1w = os.path.join(out_dir, 'sub-' + sub + '_run-01_T1w.nii.gz')

os.system('cp ' + T1map + ' ' + out_dir + 'sub-' + sub + '_run-01_T1map.nii.gz')
T1map = os.path.join(out_dir, 'sub-' + sub + '_run-01_T1map.nii.gz')

T1w_slab = os.path.join(out_dir, 'sub-' + sub + '_run-02_T1w.nii.gz')

############################################################################
# 4. Inhomogeneity correction and skull stripping
# ----------------
# Here, we perform inhomogeneity correction and skull stripping using SPM's
# segment routine. Based on the WM, GM, and CSF segmentation a brainmask is
# created and applied to the inhomogeneity corrected volumes.
print('')
print('*****************************************************')
print('* Inhomogeneity correction and skull stripping of whole brain data')
print('* Started at: ' + strftime("%Y-%m-%d %H:%M:%S", localtime()))
print('*****************************************************')

# Check for files and hires processing. Start MATLAB and run SPM with custom parameters.
if hires == True:
	checkFile = 'sub-' + sub + '_run-01_T1map_resampled_biasCorrected_masked.nii.gz'
else:
	checkFile = 'sub-' + sub + '_run-01_T1map_biasCorrected_masked.nii.gz'

if os.path.isfile(os.path.join(out_dir, checkFile)) and reprocess != True:
	print('File exists already. Skipping process.')
else:
	os.system('matlab -nosplash -nodisplay -r \"preproc_sensemap(\'' + T1w + '\'); exit;\"')

	os.system('cp ' + T1map + ' ' + out_dir + 'sub-' + sub + '_run-01_T1map.nii.gz')

	os.system('matlab -nosplash -nodisplay -r \"preproc_sensemap(\'' + T1map + '\'); exit;\"')

if hires == True:
	print('')
	print('*****************************************************')
	print('* Inhomogeneity correction and skull stripping of slab data')
	print('* Started at: ' + strftime("%Y-%m-%d %H:%M:%S", localtime()))
	print('*****************************************************')
	# Copy data and update file name
	os.system('cp ' + T1map_slab + ' ' + out_dir + 'sub-' + sub + '_run-02_T1map.nii.gz')
	T1map_slab = os.path.join(out_dir, 'sub-' + sub + '_run-02_T1map.nii.gz')

	if os.path.isfile(os.path.join(out_dir, 'sub-' + sub + '_run-02_T1map_biasCorrected.nii.gz')) and reprocess != True:
		print('File exists already. Skipping process.')
	else:
		os.system('matlab -nosplash -nodisplay -r \"preproc_sensemap(\'' + T1w_slab + '\'); exit;\"')
		os.system('matlab -nosplash -nodisplay -r \"preproc_sensemap(\'' + T1map_slab + '\'); exit;\"')

	# Update file names
	mask = os.path.join(out_dir, 'sub-' + sub + '_run-01_T1map_resampled_brainmask.nii.gz')
	T1w_biasCorrected = os.path.join(out_dir, 'sub-' + sub + '_run-01_T1w_resampled_biasCorrected.nii.gz')
	T1w_biasCorrected_masked = os.path.join(out_dir, 'sub-' + sub + '_run-01_T1w_resampled_biasCorrected_masked.nii.gz')
	T1map_biasCorrected = os.path.join(out_dir, 'sub-' + sub + '_run-01_T1map_resampled_biasCorrected.nii.gz')
	T1map_biasCorrected_masked = os.path.join(out_dir, 'sub-' + sub + '_run-01_T1map_resampled_biasCorrected_masked.nii.gz')

	T1map_slab_biasCorrected = os.path.join(out_dir, 'sub-' + sub + '_run-02_T1map_biasCorrected.nii.gz')
	T1w_slab_biasCorrected = os.path.join(out_dir, 'sub-' + sub + '_run-02_T1w_biasCorrected.nii.gz')
	T1map_slab_biasCorrected_masked = os.path.join(out_dir, 'sub-' + sub + '_run-02_T1map_biasCorrected_masked.nii.gz')
else:
	# Mask uncorrected T1map if it does not exist already
	print('')
	print('*****************************************************')
	print('* Masking T1map')
	print('* Started at: ' + strftime("%Y-%m-%d %H:%M:%S", localtime()))
	print('*****************************************************')
	# Update file names	
	T1w_biasCorrected = os.path.join(out_dir, 'sub-' + sub + '_run-01_T1w_biasCorrected.nii.gz')
	T1w_biasCorrected_masked = os.path.join(out_dir, 'sub-' + sub + '_run-01_T1w_biasCorrected_masked.nii.gz')
	T1map_biasCorrected = os.path.join(out_dir, 'sub-' + sub + '_run-01_T1map_biasCorrected.nii.gz')
	T1map_biasCorrected_masked = os.path.join(out_dir, 'sub-' + sub + '_run-01_T1map_biasCorrected_masked.nii.gz')
	mask = os.path.join(out_dir, 'sub-' + sub + '_run-01_T1map_brainmask.nii.gz')
	
	# Check if volume exists
	if os.path.isfile(os.path.join(out_dir, 'sub-' + sub + '_run-01_T1map_masked.nii.gz')) and reprocess != True:
		print('File exists already. Skipping process.')
	else:
	# Mask image otherwise
		maskedImage = ants.mask_image(ants.image_read(T1map), ants.image_read(mask))
		ants.image_write(maskedImage, os.path.join(out_dir, 'sub-' + sub + '_run-01_T1map_masked.nii.gz'))
		print('Done.')

	# Update file names
	T1map_masked = os.path.join(out_dir, 'sub-' + sub + '_run-01_T1map_masked.nii.gz')

############################################################################
# 6.1 Register additional data to (upsampled) T1map
# -------------------
if reprocess_map_data:
	reprocess = True

if map_file_onto_surface:
	if hires == True:
		print('')
		print('*****************************************************')
		print('* Register additional data to upsampled T1map')
		print('* Started at: ' + strftime("%Y-%m-%d %H:%M:%S", localtime()))
		print('*****************************************************')
		reg1 = 'sub-' + sub + '_' + map_data_output + '_registered_to_sub-' + sub + '_run-01_T1map_resampled_biasCorrected'
	else:
		print('')
		print('*****************************************************')
		print('* Register additional data to T1map')
		print('* Started at: ' + strftime("%Y-%m-%d %H:%M:%S", localtime()))
		print('*****************************************************')
		reg1 = 'sub-' + sub + '_' + map_data_output + '_registered_to_sub-' + sub + '_run-01_T1map_biasCorrected'

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
			moving = ants.image_read(T1map_biasCorrected_masked),
			type_of_transform = 'SyNRA',
			aff_random_sampling_rate = 1,
			grad_step = 0.1,
			reg_iterations = (200, 100, 50 , 30),
			verbose = True,
			# outprefix = out_dir + reg1 + '_',
			)

			# Apply transformation
			warpedImage = ants.apply_transforms(
			fixed = ants.image_read(T1map_biasCorrected),
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
# 6.2. Apply transformation to data to be mapped on the surface
# -------------------
# Define file name for output of transformation.
if map_file_onto_surface and map_transform_file_onto_surface:
	if isinstance(transform_data_output, list):
		reg2 = []
		for name in transform_data_output:
			if hires == True:
				reg2.append('sub-' + sub + '_' + name + '_registered_to_' + sub + '_run-01_T1map_resampled_biasCorrected.nii.gz')
			else:
				reg2.append('sub-' + sub + '_' + name + '_registered_to_' + sub + '_run-01_T1map_biasCorrected.nii.gz')

	else:
		if hires == True:
			reg2 = 'sub-' + sub + '_' + transform_data_output + '_registered_to_' + sub + '_run-01_T1map_resampled_biasCorrected.nii.gz'
		else:
			reg2 = 'sub-' + sub + '_' + transform_data_output + '_registered_to_' + sub + '_run-01_T1map_biasCorrected.nii.gz'

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
						fixed = ants.image_read(T1map_biasCorrected),
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
					fixed = ants.image_read(T1map_biasCorrected),
					moving = ants.image_read(transform_data),
				        transformlist = registeredImage['invtransforms'],
					interpolator = interpolation_method,
					#imagetype = 3,
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

if reprocess_map_data:
	reprocess = False

print('')
print('*****************************************************')
print('* Processing finished successfully at: ' + strftime("%Y-%m-%d %H:%M:%S", localtime()))
print('*****************************************************')