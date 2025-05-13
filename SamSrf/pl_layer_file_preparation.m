function pl_layer_file_preparation(RootPath, DataPath, SubjDir, Subject, Subj, Condition, Order, Run, SPM, PRF, Dev, Sub, ROI)
% Input Order: RootPath, SubjPath, SubjDir, FuncPath, Subject, Condition, RunOrder, Session, Order, Run, SPM, PRF

% Prepare functional data, aperture and freesurfer folders for modelling

% .........................................................................
% Inputs
% .........................................................................
% RootPath: Modelling path [cell]
% SubjPath: Freesurfer subject path [cell]
% SubjDir: Subject directory [cell]
% FuncPath: Subject functional data directory [cell]
% Subject: Subject [cell]
% Condition: The stimulation condition (D2 and/or D3 and/or D2+D3) [cell]
% RunOrder: The stimulation order ('K1_K4', 'K4_K1') [cell]
% Session: The stimulation session ('sess1', 'sess2') [cell]
% Order: The stimulation order ('forward' and/or 'backward') [cell]
% Run: The run number ('01' and/or '02') [cell]
% SPM: FldSPM [cell]
% PRF: FldPRF [cell]
% .........................................................................
% Outputs
%..........................................................................

% Written by P.Liu
% Optimized by P.Liu
% Email: peng.liu@med.ovgu.de
% Last Updated 11 Oct 2022
%% ........................................................................Function

DirSPM = fullfile(SubjDir,[SPM Condition]);
mkdir (DirSPM);

DirPRF = fullfile(SubjDir,[PRF Condition]);
mkdir (DirPRF);

copyfile(fullfile(RootPath,'aps_verbars_vec.mat'), DirPRF);

for i_order = 1:size(Order, 2)
    
    CurrOrder = Order{i_order};
    
    for i_run = 1:size(Run, 2)
        
        CurrRun = Run{i_run};
        
        SourceDir = fullfile(DataPath, Subject, Dev, [Sub Subj]);
        
        FunImg = ([Sub Subj '_' Condition '_' CurrOrder CurrRun '_registered_to_' Subj '_run-01_T1map.nii.gz']);
        
        FileName = [CurrOrder CurrRun '.nii.gz'];
        
        copyfile(fullfile(SourceDir,FunImg),fullfile(DirSPM,FileName));
        
        copyfile(fullfile(SourceDir,ROI),fullfile(DirSPM,ROI));
        
    end
        
end

