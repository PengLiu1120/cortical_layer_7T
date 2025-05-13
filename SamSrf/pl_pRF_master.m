%% Population receptive field (pRF) modelling
% This is the master script to run pRF modelling using SamSrf
%% Functions called in this script:
% pl_file_preparation
% ss_samsrf_mgh2srf
% ss_samsrf_makesemoroi
% pl_pRF_concatenate
% pl_pRF_average
% pl_copyfile
% pl_Gaussian_Tuning_Curve_Vertical
% pl_copylabelregister
% pl_result_extraction
%% ........................................................................
% Written by P.Liu
% Optimized by P.Liu
% Email: peng.liu@med.ovgu.de
% Last updated 11 Apr 2022 by P.Liu
%% ........................................................................Tidy up
clear all
close all
clc

%% ........................................................................Get path to SamSrf
addpath(genpath('/Users/pliu/Documents/Programmes/samsrf_v9.4'));
%% ........................................................................Set defaults
% .........................................................................Specify modelling path
RootPath = '/Volumes/IKND/AG_Kuehn/Peng/LayerPRF/pRFmodel';

% .........................................................................Specify freesurfer subject path
SubjPath = '/Users/pliu/Documents/LayerPRF/subjects';
SessPath = '/Users/pliu/Documents/LayerPRF/sessions';

% .........................................................................Specify functional data path
DataPath = '/Volumes/IKND/AG_Kuehn/Peng/LayerPRF/LayerPRF';
Functional = 'Functional';
RunOrder = {'K1_K4' 'K4_K1'};
Session = {'sess1' 'sess2'};

% .........................................................................Specify subjects
Subjects = {'ajz367' 'bkn792' 'bmg520' 'cxc075' 'czg996' 'ejk164' 'frj712' 'ggp057' 'gph998' 'gxo876' 'hby152' 'ijt563' 'iwq192' 'kdy341' 'llh150' 'lpr469' 'nhm378' 'oms448' 'qet940' 'qxo538' 'sst050' 'unk742'};
% 'ajz367' 'bkn792' 'bmg520' 'cxc075' 'czg996' 'ejk164' 'frj712' 'ggp057' 'gph998' 'gxo876' 'hby152' 'hji464' 'ijt563' 'iwq192' 'kdy341' 'llh150' 'lpr469' 'nhm378' 'oms448' 'qet940' 'qxo538' 'sst050' 'unk742'

% .........................................................................Specify register.dat path and labels
Label = 'label';
Image = 'image';
Fourier = '_Fourier';
FuncSession = '_backward_1';

% .........................................................................Specify condtions, D2, D3 and D2+D3
Conditions = {'D2' 'D3' 'D2+D3'};
% 'D2' 'D3' 'D2+D3'

% .........................................................................Specify hemisphere
Hemis = {'lh' 'rh'};

% .........................................................................Specify localisers, hand and 3b
Loc = {'3b_hand'};
% '3b_hand' '3b'

% .........................................................................Specify parameters to convert functional data from mgh format (previously converted from nifti format by Freesurfer) to srf format.
Combi = 'runwise';
FldSrf = 'surf';
FldSPM = 'spm_';
FldPRF = 'pRF_';
FuncPattern = '*h_*.mgh';
Norm = true;
Avg = false;
NoiseCeiling = false;

% .........................................................................Specify parameters to concatenate and average files generated from mgh2srf
Order = {'forward' 'backward'};
RunConc = {'01' '02'};
RunAvg = {'forward01backward01' 'forward02backward02'};

% .........................................................................Specify Model Fit parameters
ModelFiles = 'Average';
ROIFiles = 'sensorimotor';

% .........................................................................Switch for functions
% .........................................................................Step One
switch_prepare = false;
switch_copyregister = false;

% .........................................................................Generate MGH files
% .........................................................................Shell script
        
% .........................................................................Step Two
switch_mgh2srf = false;
switch_createROI = false;
switch_conc_avg = false;
switch_copy = false;
switch_modelfit = false;
% .........................................................................Step Three
switch_label = false;
switch_results = true;

%% ........................................................................Preprocessing pipeline
for i_sub=1:size(Subjects, 2)
    
    CurrSubj = Subjects{i_sub};
    
    % .....................................................................Make subject directory
    CurrSubjPath = fullfile(RootPath,CurrSubj);
    mkdir (CurrSubjPath);
    
    % .....................................................................Make subject surf directory
    CurrPathSrf = fullfile(CurrSubjPath,FldSrf);
    
    for i_cond = 1:size(Conditions,2)
        
        CurrCond = Conditions{i_cond};
        
        % .................................................................Make subject functional data directory
        CurrDataPath = fullfile(DataPath,CurrSubj,Functional);
        CurrFuncPath = fullfile(CurrDataPath, CurrCond);
        
        % .................................................................Make directory and copy files NIFTI
        
        if switch_prepare
            
            pl_file_preparation(RootPath, SubjPath, CurrSubjPath, CurrFuncPath, CurrSubj, CurrCond, RunOrder, Session, Order, RunConc, FldSPM, FldPRF);
            
        end
        
        CurrPathSPM = fullfile(CurrSubjPath,[FldSPM CurrCond]);
        CurrPathPRF = fullfile(CurrSubjPath,[FldPRF CurrCond]);
        
        % .................................................................Copy register.dat file
        
        if switch_copyregister
            
            pl_copyregister(SessPath, CurrPathSPM, CurrSubj, CurrCond, Image, Fourier, FuncSession)
            
        end
        
        % .................................................................mgh2srf
        
        if switch_mgh2srf
            
            cd(CurrPathSPM)
            
            ss_samsrf_mgh2srf(Hemis, Combi, CurrPathSrf, FuncPattern, Norm, Avg, NoiseCeiling);
            
        end
        
        % .................................................................CreateROI
        
        if switch_createROI
            
            cd(CurrPathPRF)
            
            ss_samsrf_makesemoroi(CurrPathSrf);
            
        end
        
        % .................................................................Concatenate and average
        
        if switch_conc_avg
            
            cd(CurrPathSPM)
            
            pl_pRF_concatenate(Hemis, Order, RunConc)
            pl_pRF_average(Hemis, RunAvg)
            
        end
        
        % .................................................................Copy over ModelFit files
        
        if switch_copy
            
            pl_copyfile(Hemis, ModelFiles, CurrPathSPM, CurrPathPRF)
            
        end
        
        % .................................................................Gaussian Tuning Curve Vertical
        
        if switch_modelfit
            
            cd(CurrPathPRF)
            
            pl_Gaussian_Tuning_Curve_Vertical(Hemis, ModelFiles, ROIFiles)
            
        end
        
        % .................................................................Copy 3b hand & localiser label
        
        if switch_label
            
            for i_Loc = 1:size(Loc,2)
                
                CurrLoc = Loc{i_Loc};
                
                pl_copylabel(SubjPath, CurrPathPRF, CurrSubj, Label, CurrLoc)
                
            end
            
        end
        
        % .................................................................Results extraction
        
        if switch_results
            
            cd (CurrPathPRF)
            
            for i_Loc = 1:size(Loc,2)
                
                CurrLoc = Loc{i_Loc};
            
                % pl_result_extraction(Hemis, CurrCond)
                pl_result_extraction(CurrCond, CurrLoc)
                
            end
            
        end
        
    end
    
end

%% ........................................................................Remove SamSrf path
rmpath(genpath('/Users/pliu/Documents/Programmes/samsrf_v9.4'));