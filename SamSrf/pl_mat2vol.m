function pl_mat2vol(SubjDir, PRF, Condition)
% Input order: 

% Converts SamSrf .mat files to Nifti

% .........................................................................
% Inputs
% .........................................................................

% .........................................................................
% Outputs
% .........................................................................
% Files can be converted including:
% Raw data converted from Nifti to .mat
% Concatnated or averaged .mat file
% .mat file after fitting
% .........................................................................
% Note: Converted volumetric data will be in 3D space
% .........................................................................
% Written by P.Liu
% Optimized by P.Liu
% Email: peng.liu@med.ovgu.de
% Last updated 17 Oct 2022 by P.Liu
%% ........................................................................Function

DirPRF = fullfile(SubjDir,[PRF Condition]);
cd(DirPRF)

SrfName = 'vol_pTC_ver';

samsrf_mat2vol(SrfName)

end