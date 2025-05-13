%% ........................................................................Get path to SamSrf
addpath(genpath('/Users/pliu/Documents/Programmes/samsrf_v9.4'));

label = 'lh.3b_hand';
funimg = 'sub-frj_task-phase_bold_D2+D3_mean';
strimg = 'T1';
hemi = 'lh';

samsrf_label2nii(label, funimg, strimg, hemi)

%% ........................................................................Remove SamSrf path
rmpath(genpath('/Users/pliu/Documents/Programmes/samsrf_v9.4'));