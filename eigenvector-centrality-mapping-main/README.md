# Eigenvector-Centrality-Mapping

## Requirements
These bash scripts depend on the following software:
- https://github.com/lipsia-fmri/lipsia
- FSL

Make sure to install them first. 

input data: 4D nifti times series (corrected for slice timing, motion and physiological noise)
see ECM_7T_instructions.pdf for pre-processing steps (scripts included in https://gitlab.com/estherkuehnneuroscience/fmri_preprocessing and https://gitlab.com/estherkuehnneuroscience/physiodata_preprocessing)

## Calculate maps

1. run bash script _fsl_avg_masking_updated_ to calculate time series average image (fsl5.0-fslmaths data.nii -Tmean data_Tmean.nii) and
to extract a brain mask (fsl5.0-bet2 data_Tmean.nii data_Tmean_brain -f 0.1 -g 0 -m)
1. convert mat file containg moco params (*MoCoParam*.mat) to text, for example using matlab:
    - load('moco.mat')
    - dlmwrite('moco.txt', myFile, 'delimiter','\t')
1. run bash script _vnifti-vpreprocess-vecm-final_

### vnifti-vpreprocess-vecm-final - what the script does:
calls lipsia functions to perform calculations on resting state fMRI data:
- regressing out motion parameters
- calculating ecm with rlc metric developed by Gaby Lohmann
- calculating ecm with add metric
- showing results in vini viewer

input data: fMRI time series nii (data.nii or adata.nii)
output data: covariates (motion), residuals, highpass filtered and smoothed data, ecm maps for two metrics (rlc, add)

## Registration and surface mapping of Eigenvector Centrality data:
1. Register functional statistical maps to the qT1 image or other contrast of interest:
    - go to ITK Snap (v3.8.0) > load prepared high-res qT1 slab image as reference image
    - load 3D EPI (i.e. 4D averaged across timepoints) as additional image
    - open Registration Tool
    - Functional images were registered to the qT1 slab images using the automated registration tool (non-rigid, 9 degrees of freedom)  with manual refinement (prioritizing alignment in M1 and S1) where necessary (mismatch < 1voxel).
    - resulting ITK registration matrices (stored as txt) were applied to statistical maps (i.e. pRF maps, tmaps) using ANTs v2.1.0 by  running the following bash script: ecm_map_preprocess
1. Map registered data to cortical surfaces:
    - open MIPAV and run pipeline Mapping_function_ecm.LayoutXML 
    - https://gitlab.com/estherkuehnneuroscience/mipav/-/blob/master/pipelines/3D_Architecture_S1_Aging/Mapping_function_ecm.LayoutXML?ref_type=heads
