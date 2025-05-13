function pl_vol2mat(SubjDir, SPM, Condition, Order, Run, ROImg, ROI)
% Input order: SubjDir, SPM, Condition, Order, Run

% Convert volumetric functional data (3D/4D) to SamSrf .mat file 

% .........................................................................
% Inputs
% .........................................................................
% SubjDir: Subject directory [cell]
% SPMPath: Path to spm_ folder [cell]
% Condition: The stimulation condition (D2 and/or D3 and/or D2+D3) [cell]
% Order: The stimulation order (forward and backward) [cell]
% Run: The stimulation run (01 and 02) [cell]
% .........................................................................
% Outputs
% .........................................................................
% Written by P.Liu
% Optimized by P.Liu
% Email: peng.liu@med.ovgu.de
% Last updated 11 Oct 2022 by P.Liu
%% ........................................................................Function

DirSPM = fullfile(SubjDir,[SPM Condition]);
cd(DirSPM)

% gunzip(ROImg);
% delete(ROImg);

for i_order = 1:size(Order, 2)
    
    CurrOrder = Order{i_order};
    
    for i_run = 1:size(Run, 2)
        
        CurrRun = Run{i_run};
        
        % FuncImg = [CurrOrder CurrRun '.nii.gz'];
        % gunzip(FuncImg);
        % delete(FuncImg);
        
        FunImg = [CurrOrder CurrRun];
        
        samsrf_vol2mat(FunImg, ROI, false);
        
    end
    
end

end