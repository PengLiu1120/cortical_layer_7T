%% Localiser data copy
% This script copies functional localiser GLM results from preprocessing
% folder to Nighres surface mapping pipeline with BID format
%% ........................................................................
% Written by P.Liu
% Optimized by P.Liu
% Email: peng.liu@med.ovgu.de
% Last updated 21 Mar 2023 by P.Liu
%% ........................................................................Tidy up
clear all
close all
clc

%% ........................................................................Set defaults
%..........................................................................Specify modelling path
RootPath = '/Volumes/IKND/AG_Kuehn/Peng/LayerPRF/LayerPRF';
DataPath = '/Volumes/IKND/AG_Kuehn/Peng/LayerPRF/SurfaceRegister';

% .........................................................................Specify subjects
Subjects = {'ajz367' 'bkn792' 'bmg520' 'bmz426' 'cxc075' 'czg996' 'ejk164' 'frj712' 'ggp057' 'gph998' 'gxo876' 'hby152' 'ijt563' 'iwq192' 'kdy341' 'llh150' 'lpr469' 'nhm378' 'oms448' 'qet940' 'qxo538' 'sst050' 'ugn780' 'unk742'};
Subj = {'ajz' 'bkn' 'bmg' 'bmz' 'cxc' 'czg' 'ejk' 'frj' 'ggp' 'gph' 'gxo' 'hby' 'ijt' 'iwq' 'kdy' 'llh' 'lpr' 'nhm' 'oms' 'qet' 'qxo' 'sst' 'ugn' 'unk'};

% .........................................................................Specify condtions, D2, D3 and D2+D3
Conditions = {'D2' 'D3' 'D2+D3'};

% .........................................................................Specify fields
FldDev = 'derivatives';
FldSub = 'sub-';
FldFun = 'Functional';
FldLoc = 'Localiser';
FldRes = 'Results';

%% ........................................................................Copy file
for i_sub=1:size(Subjects, 2)
    
    CurrSubj = Subjects{i_sub};
    CurrSub = Subj{i_sub};
    
    SourceDir_1 = fullfile(RootPath, CurrSubj, FldRes, FldLoc);
    SourceDir_2 = fullfile(RootPath, CurrSubj, FldFun, FldLoc);
    TargetDir_1 = fullfile(DataPath, CurrSubj, FldDev, [FldSub CurrSub], FldLoc);
    mkdir (TargetDir_1);
    TargetDir_2 = fullfile(DataPath, CurrSubj, FldDev, [FldSub CurrSub]);
    
    for i_cond = 1:size(Conditions,2)
        
        CurrCond = Conditions{i_cond};
        
        FunImg = [CurrCond '.nii'];
        LocImg = [CurrCond '_Loc.nii' ];
        
        copyfile(fullfile(SourceDir_1,FunImg),fullfile(TargetDir_1, LocImg));
        
    end
    
    meanImg = 'adata_mean_average.nii';
    LocmeanImg = 'localiser_mean.nii';
    
    copyfile(fullfile(SourceDir_2,meanImg),fullfile(TargetDir_2, LocmeanImg));
    
end