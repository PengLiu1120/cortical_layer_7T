function pl_layer_copyfile(ModelFiles, PathSPM, PathPRF)
% Input order: Hemis, ModelFiles, PathSPM, PathPRF

% Copy Averaged file from spm folder to prf folder

% .........................................................................
% Inputs
% .........................................................................
% Hemis: Hemispheres ('lh' and/or 'rh') [cell]
% ModelFiles: Cell array with SamSrf data files (without extension) [cell]
% PathSPM: Path to spm folder [cell]
% PathPRF: Path to prf folder [cell]
% .........................................................................
% Outputs
% .........................................................................
% Written by P.Liu
% Optimized by P.Liu
% Email: peng.liu@med.ovgu.de
% Last updated 11 Oct 2022 by P.Liu
%% ........................................................................Function
FileName = ['vol_' ModelFiles '.mat'];

copyfile(fullfile(PathSPM, FileName), PathPRF);
end