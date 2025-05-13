function pl_file_preparation(RootPath, SubjPath, SubjDir, FuncPath, Subject, Condition, RunOrder, Session, Order, Run, SPM, PRF)
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
% Last Updated 12 May 2022
%% ........................................................................Function

DirSPM = fullfile(SubjDir,[SPM Condition]);
mkdir (DirSPM);

DirPRF = fullfile(SubjDir,[PRF Condition]);
mkdir (DirPRF);

DirSubj = fullfile(SubjPath,Subject);

copyfile(DirSubj, SubjDir);

copyfile(fullfile(RootPath,'aps_verbars_vec.mat'), DirPRF);

for i_order = 1:size(Order, 2)
    
    CurrOrder = Order{i_order};
    CurrRunOrder = RunOrder{i_order};
    
    for i_session = 1:size(Session,2)
        
        CurrSession = Session{i_session};
        CurrRun = Run{i_session};
        
        SourceDir = fullfile(FuncPath, CurrRunOrder, CurrSession);
        
        FileName = [CurrOrder CurrRun '.nii'];
        
        copyfile(fullfile(SourceDir,'adata.nii'),fullfile(DirSPM,FileName));
        
    end
    
end

end