mkdir /Users/pliu/Documents/LayerPRF/PRFmodelling/frj712_hires

mkdir /Users/pliu/Documents/LayerPRF/PRFmodelling/frj712_hires/pRF_d2+d3
mkdir /Users/pliu/Documents/LayerPRF/PRFmodelling/frj712_hires/spm_d2+d3
mkdir /Users/pliu/Documents/LayerPRF/PRFmodelling/frj712_hires/pRF_d2
mkdir /Users/pliu/Documents/LayerPRF/PRFmodelling/frj712_hires/spm_d2
mkdir /Users/pliu/Documents/LayerPRF/PRFmodelling/frj712_hires/pRF_d3
mkdir /Users/pliu/Documents/LayerPRF/PRFmodelling/frj712_hires/spm_d3

cp -R /Users/pliu/Documents/LayerPRF/subjects/frj712_hires/. /Users/pliu/Documents/LayerPRF/PRFmodelling/frj712_hires

cp /Users/pliu/Documents/Programmes/SamSrf/aps_verbars.mat /Users/pliu/Documents/LayerPRF/PRFmodelling/frj712_hires/pRF_d2+d3
cp /Users/pliu/Documents/Programmes/SamSrf/aps_verbars.mat /Users/pliu/Documents/LayerPRF/PRFmodelling/frj712_hires/pRF_d2
cp /Users/pliu/Documents/Programmes/SamSrf/aps_verbars.mat /Users/pliu/Documents/LayerPRF/PRFmodelling/frj712_hires/pRF_d3

cp /Users/pliu/Documents/LayerPRF/LayerPRF/frj712_hires/Functional/D2+D3/K1_K4/sess1/adata.nii /Users/pliu/Documents/LayerPRF/PRFmodelling/frj712_hires/spm_d2+d3/forward01.nii
cp /Users/pliu/Documents/LayerPRF/LayerPRF/frj712_hires/Functional/D2+D3/K1_K4/sess2/adata.nii /Users/pliu/Documents/LayerPRF/PRFmodelling/frj712_hires/spm_d2+d3/forward02.nii
cp /Users/pliu/Documents/LayerPRF/LayerPRF/frj712_hires/Functional/D2+D3/K4_K1/sess1/adata.nii /Users/pliu/Documents/LayerPRF/PRFmodelling/frj712_hires/spm_d2+d3/backward01.nii
cp /Users/pliu/Documents/LayerPRF/LayerPRF/frj712_hires/Functional/D2+D3/K4_K1/sess2/adata.nii /Users/pliu/Documents/LayerPRF/PRFmodelling/frj712_hires/spm_d2+d3/backward02.nii

cp /Users/pliu/Documents/LayerPRF/LayerPRF/frj712_hires/Functional/D2/K1_K4/sess1/adata.nii /Users/pliu/Documents/LayerPRF/PRFmodelling/frj712_hires/spm_d2/forward01.nii
cp /Users/pliu/Documents/LayerPRF/LayerPRF/frj712_hires/Functional/D2/K1_K4/sess2/adata.nii /Users/pliu/Documents/LayerPRF/PRFmodelling/frj712_hires/spm_d2/forward02.nii
cp /Users/pliu/Documents/LayerPRF/LayerPRF/frj712_hires/Functional/D2/K4_K1/sess1/adata.nii /Users/pliu/Documents/LayerPRF/PRFmodelling/frj712_hires/spm_d2/backward01.nii
cp /Users/pliu/Documents/LayerPRF/LayerPRF/frj712_hires/Functional/D2/K4_K1/sess2/adata.nii /Users/pliu/Documents/LayerPRF/PRFmodelling/frj712_hires/spm_d2/backward02.nii

cp /Users/pliu/Documents/LayerPRF/LayerPRF/frj712_hires/Functional/D3/K1_K4/sess1/adata.nii /Users/pliu/Documents/LayerPRF/PRFmodelling/frj712_hires/spm_d3/forward01.nii
cp /Users/pliu/Documents/LayerPRF/LayerPRF/frj712_hires/Functional/D3/K1_K4/sess2/adata.nii /Users/pliu/Documents/LayerPRF/PRFmodelling/frj712_hires/spm_d3/forward02.nii
cp /Users/pliu/Documents/LayerPRF/LayerPRF/frj712_hires/Functional/D3/K4_K1/sess1/adata.nii /Users/pliu/Documents/LayerPRF/PRFmodelling/frj712_hires/spm_d3/backward01.nii
cp /Users/pliu/Documents/LayerPRF/LayerPRF/frj712_hires/Functional/D3/K4_K1/sess2/adata.nii /Users/pliu/Documents/LayerPRF/PRFmodelling/frj712_hires/spm_d3/backward02.nii

cp /Users/pliu/Documents/LayerPRF/sessions/image/frj712_hires_Fourier/image/D2+D3_backward_1/register.dat /Users/pliu/Documents/LayerPRF/PRFmodelling/frj712_hires/spm_d2+d3/register.dat
cp /Users/pliu/Documents/LayerPRF/sessions/image/frj712_hires_Fourier/image/D2_backward_1/register.dat /Users/pliu/Documents/LayerPRF/PRFmodelling/frj712_hires/spm_d2/register.dat
cp /Users/pliu/Documents/LayerPRF/sessions/image/frj712_hires_Fourier/image/D3_backward_1/register.dat /Users/pliu/Documents/LayerPRF/PRFmodelling/frj712_hires/spm_d3/register.dat

cp /Users/pliu/Documents/LayerPRF/subjects/frj712_hires/label/lh.3b_hand_loc_D2.label /Users/pliu/Documents/LayerPRF/PRFmodelling/frj712_hires/pRF_d2/lh.3b_hand_loc_d2.label
cp /Users/pliu/Documents/LayerPRF/subjects/frj712_hires/label/lh.3b_hand_loc_D3.label /Users/pliu/Documents/LayerPRF/PRFmodelling/frj712_hires/pRF_d3/lh.3b_hand_loc_d3.label
cp /Users/pliu/Documents/LayerPRF/subjects/frj712_hires/label/lh.3b_hand_loc_D2+D3.label /Users/pliu/Documents/LayerPRF/PRFmodelling/frj712_hires/pRF_d2+d3/lh.3b_hand_loc_d2+d3.label