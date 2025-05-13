function pl_copyregister(SessPath, SPMPath, Subject, Condition, Image, Fourier, FuncSession)
% function pl_copylabelregister(Hemis, SubjPath, SessPath, SPMPath, PRFPath, Subject, Label, Condition, Image, Fourier, FuncSession)
% Input order: Hemis, ModelFiles, PathSPM, PathPRF

% Copy labels used for data extraction
% Labels are 3b hand + localiser overlap

% Only left hemisphere

% .........................................................................
% Inputs
% .........................................................................
% Hemis: Hemispheres ('lh' and/or 'rh') [cell]
% SubjPath: Freesurfer subject path [cell]
% SessPath: Freesurfer session path [cell]
% SPMPath: Path to spm_ folder [cell]
% PRFPath: Path to prf_ folder [cell]
% Subject: Subject [cell]
% Label: label [cell]
% Condition: The stimulation condition (D2 and/or D3 and/or D2+D3) [cell]
% Image: 'image' [cell]
% Fourier: '_Fourier' [cell]
% FuncSession: '_backward_1' [cell]
% .........................................................................
% Outputs
% .........................................................................
% Written by P.Liu
% Optimized by P.Liu
% Email: peng.liu@med.ovgu.de
% Last updated 25 Mar 2022 by P.Liu
%% ........................................................................Function

% for i_hemi = 1:size(Hemis, 2)
    
    % CurrHemi = Hemis{i_hemi};
    
    % .........................................................................copy register.dat
    RegisterDir = fullfile(SessPath,Image,[Subject Fourier], Image, [Condition FuncSession]);
    
    copyfile(fullfile(RegisterDir,'register.dat'), SPMPath);
% end

end