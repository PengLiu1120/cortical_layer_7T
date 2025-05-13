%% ........................................................................Tidy up
clear all
close all
clc

%% ........................................................................Set defaults
%..........................................................................Specify modelling path
RootPath = '/Volumes/IKND/AG_Kuehn/Peng/LayerPRF/pRFmodel';
DataPath = '/Volumes/IKND/AG_Kuehn/Peng/LayerPRF/SurfaceRegister';

% .........................................................................Specify subjects
Subjects = {'ajz367' 'bkn792' 'bmg520' 'bmz426' 'cxc075' 'czg996' 'ejk164' 'frj712' 'ggp057' 'gph998' 'gxo876' 'hby152' 'ijt563' 'iwq192' 'kdy341' 'llh150' 'lpr469' 'nhm378' 'oms448' 'qet940' 'qxo538' 'sst050' 'ugn780' 'unk742'};
Subj = {'ajz' 'bkn' 'bmg' 'bmz' 'cxc' 'czg' 'ejk' 'frj' 'ggp' 'gph' 'gxo' 'hby' 'ijt' 'iwq' 'kdy' 'llh' 'lpr' 'nhm' 'oms' 'qet' 'qxo' 'sst' 'ugn' 'unk'};

% .........................................................................Specify fields
FldDev = 'derivatives';
FldSub = 'sub-';
FldSPM = 'spm_';
FldPRF = 'pRF_';

% .........................................................................Specify condtions, D2, D3 and D2+D3
Conditions = {'D2' 'D3' 'D2+D3'};

% .........................................................................Specify image name
Order = {'forward' 'backward'};
RunCon = {'01' '02'};

%% ........................................................................Copy file
for i_sub=1:size(Subjects, 2)
    
    CurrSubj = Subjects{i_sub};
    CurrSub = Subj{i_sub};
    
    for i_cond = 1:size(Conditions,2)
        
        CurrCond = Conditions{i_cond};
        
        for i_order = 1:size(Order, 2)
            
            CurrOrder = Order{i_order};
            
            for i_run = 1:size(RunCon, 2)
                
                CurrRun = RunCon{i_run};
                
                SourceDir = fullfile(RootPath, CurrSubj, [FldSPM CurrCond]);
                TargetDir = fullfile(DataPath, CurrSubj, FldDev, [FldSub CurrSub], [FldPRF CurrCond]);
                
                FunImg = [CurrOrder CurrRun '.nii'];
                
                copyfile(fullfile(SourceDir,FunImg),TargetDir);
                
            end
            
        end
        
    end
    
end