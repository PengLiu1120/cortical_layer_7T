#!/bin/bash
# pl_swapdim_afni
#
#-----------------------------------------------------------------------------
# Inputs
#-----------------------------------------------------------------------------
# $1 - Subject IDs
#-----------------------------------------------------------------------------
# Outputs
#-----------------------------------------------------------------------------
#  
# ----------------------------------------------------------------------------
# Usage:
#
#-----------------------------------------------------------------------------
# Note: 
#-----------------------------------------------------------------------------
# 
# 31/05/2022: Generated (PL)
# 31/05/2022: Last modified (PL)

#-----------------------------------------------------------------------------

SubjIDs=$1

for i_subj in $SubjIDs
do

cd /Users/pliu/Documents/LayerPRF/LayerPRF/kdy341/Functional/D2/K1_K4/sess1
fslswapdim adata_raw.nii RL SI PA adata.nii
gunzip adata.nii.gz
fslswapdim afourier_raw.nii RL SI PA afourier.nii
gunzip afourier.nii.gz
3dcopy afourier.nii afourier+orig.BRIK

cd /Users/pliu/Documents/LayerPRF/LayerPRF/kdy341/Functional/D2/K1_K4/sess2
fslswapdim adata_raw.nii RL SI PA adata.nii
gunzip adata.nii.gz
fslswapdim afourier_raw.nii RL SI PA afourier.nii
gunzip afourier.nii.gz
3dcopy afourier.nii afourier+orig.BRIK

cd /Users/pliu/Documents/LayerPRF/LayerPRF/kdy341/Functional/D2/K4_K1/sess1
fslswapdim adata_raw.nii RL SI PA adata.nii
gunzip adata.nii.gz
fslswapdim afourier_raw.nii RL SI PA afourier.nii
gunzip afourier.nii.gz
3dcopy afourier.nii afourier+orig.BRIK

cd /Users/pliu/Documents/LayerPRF/LayerPRF/kdy341/Functional/D2/K4_K1/sess2
fslswapdim adata_raw.nii RL SI PA adata.nii
gunzip adata.nii.gz
fslswapdim afourier_raw.nii RL SI PA afourier.nii
gunzip afourier.nii.gz
3dcopy afourier.nii afourier+orig.BRIK

cd /Users/pliu/Documents/LayerPRF/LayerPRF/kdy341/Functional/D2+D3/K1_K4/sess1
fslswapdim adata_raw.nii RL SI PA adata.nii
gunzip adata.nii.gz
fslswapdim afourier_raw.nii RL SI PA afourier.nii
gunzip afourier.nii.gz
3dcopy afourier.nii afourier+orig.BRIK

cd /Users/pliu/Documents/LayerPRF/LayerPRF/kdy341/Functional/D2+D3/K1_K4/sess2
fslswapdim adata_raw.nii RL SI PA adata.nii
gunzip adata.nii.gz
fslswapdim afourier_raw.nii RL SI PA afourier.nii
gunzip afourier.nii.gz
3dcopy afourier.nii afourier+orig.BRIK

cd /Users/pliu/Documents/LayerPRF/LayerPRF/kdy341/Functional/D2+D3/K4_K1/sess1
fslswapdim adata_raw.nii RL SI PA adata.nii
gunzip adata.nii.gz
fslswapdim afourier_raw.nii RL SI PA afourier.nii
gunzip afourier.nii.gz
3dcopy afourier.nii afourier+orig.BRIK

cd /Users/pliu/Documents/LayerPRF/LayerPRF/kdy341/Functional/D2+D3/K4_K1/sess2
fslswapdim adata_raw.nii RL SI PA adata.nii
gunzip adata.nii.gz
fslswapdim afourier_raw.nii RL SI PA afourier.nii
gunzip afourier.nii.gz
3dcopy afourier.nii afourier+orig.BRIK

cd /Users/pliu/Documents/LayerPRF/LayerPRF/kdy341/Functional/D3/K1_K4/sess1
fslswapdim adata_raw.nii RL SI PA adata.nii
gunzip adata.nii.gz
fslswapdim afourier_raw.nii RL SI PA afourier.nii
gunzip afourier.nii.gz
3dcopy afourier.nii afourier+orig.BRIK

cd /Users/pliu/Documents/LayerPRF/LayerPRF/kdy341/Functional/D3/K1_K4/sess2
fslswapdim adata_raw.nii RL SI PA adata.nii
gunzip adata.nii.gz
fslswapdim afourier_raw.nii RL SI PA afourier.nii
gunzip afourier.nii.gz
3dcopy afourier.nii afourier+orig.BRIK

cd /Users/pliu/Documents/LayerPRF/LayerPRF/kdy341/Functional/D3/K4_K1/sess1
fslswapdim adata_raw.nii RL SI PA adata.nii
gunzip adata.nii.gz
fslswapdim afourier_raw.nii RL SI PA afourier.nii
gunzip afourier.nii.gz
3dcopy afourier.nii afourier+orig.BRIK

cd /Users/pliu/Documents/LayerPRF/LayerPRF/kdy341/Functional/D3/K4_K1/sess2
fslswapdim adata_raw.nii RL SI PA adata.nii
gunzip adata.nii.gz
fslswapdim afourier_raw.nii RL SI PA afourier.nii
gunzip afourier.nii.gz
3dcopy afourier.nii afourier+orig.BRIK

done
