function pl_copyfile(Hemis, ModelFiles, PathSPM, PathPRF)
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
% Last updated 11 Apr 2022 by P.Liu
%% ........................................................................Function

% .........................................................................Loop through Hemis
for i_hemi = 1:size(Hemis, 2)
    
    CurrHemi = Hemis{i_hemi};
    
    FileName = [CurrHemi '_' ModelFiles '.mat'];
    
    copyfile(fullfile(PathSPM, FileName), PathPRF);
end

end