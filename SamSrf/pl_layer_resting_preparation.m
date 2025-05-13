function pl_layer_resting_preparation(DataPath, SubjDir, Subject, Subj, Condition, SPM, PRF, Dev, Sub)
% Input Order: DataPath, SubjPath, SubjDir, Subject, Condition, SPM, PRF, Dev, Sub

% Prepare resting state images for BOLD signal correction

% .........................................................................
% Inputs
% .........................................................................
% RootPath: Modelling path [cell]
% SubjDir: Subject directory [cell]
% Subject: Subject [cell]
% Condition: The stimulation condition (D2 and/or D3 and/or D2+D3) [cell]
% SPM: FldSPM [cell]
% PRF: FldPRF [cell]
% Dev: FldDev [cell]
% Sub: Subj [cell]
% .........................................................................
% Outputs
%..........................................................................

% Written by P.Liu
% Optimized by P.Liu
% Email: peng.liu@med.ovgu.de
% Last Updated 13 Mar 2023
%% ........................................................................Function

DirSPM = fullfile(SubjDir,[SPM Condition]);
mkdir (DirSPM);

DirPRF = fullfile(SubjDir,[PRF Condition]);
mkdir (DirPRF);

SourceDir = fullfile(DataPath, Subject, Dev, [Sub Subj]);

FunImg = ([Sub Subj '_resting_state_registered_to_' Subj '_run-01_T1map.nii.gz']);

FileName = 'resting_state.nii.gz';

copyfile(fullfile(SourceDir,FunImg),fullfile(DirSPM,FileName));
                
end

