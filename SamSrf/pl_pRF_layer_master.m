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
% Last updated 17 Oct 2022 by P.Liu
%% ........................................................................Tidy up
clear all
close all
clc

%% ........................................................................Get path to SamSrf
addpath(genpath('/Users/pliu/Documents/Programmes/samsrf_v9.4'));

%% ........................................................................Set defaults
%..........................................................................Specify modelling path
RootPath = '/Volumes/IKND/AG_Kuehn/Peng/LayerPRF/LayerpRFmodel';
DataPath = '/Volumes/IKND/AG_Kuehn/Peng/LayerPRF/SurfaceRegister';

% .........................................................................Specify freesurfer subject path
SubjPath = '/Users/pliu/Documents/LayerPRF/subjects';
SessPath = '/Users/pliu/Documents/LayerPRF/sessions';

% .........................................................................Specify subjects
Subjects = {'ajz367' 'bkn792' 'bmg520' 'bmz426' 'cxc075' 'czg996' 'ejk164' 'frj712' 'ggp057' 'gph998' 'gxo876' 'hby152' 'ijt563' 'iwq192' 'kdy341' 'llh150' 'lpr469' 'nhm378' 'oms448' 'qet940' 'qxo538' 'sst050' 'unk742'};
Subj = {'ajz' 'bkn' 'bmg' 'bmz' 'cxc' 'czg' 'ejk' 'frj' 'ggp' 'gph' 'gxo' 'hby' 'ijt' 'iwq' 'kdy' 'llh' 'lpr' 'nhm' 'oms' 'qet' 'qxo' 'sst' 'unk'};

% .........................................................................Specify fields
FldDev = 'derivatives';
FldSub = 'sub-';
FldSPM = 'spm_';
FldPRF = 'pRF_';

% .........................................................................Specify condtions, D2, D3 and D2+D3
Conditions = {'D2' 'D3' 'D2+D3'};
% 'D2' 'D3' 'D2+D3'

% .........................................................................Specify image name
Order = {'forward' 'backward'};
RunCon = {'01' '02'};
RunAvg = {'forward01backward01' 'forward02backward02'};

% .........................................................................Specify register.dat path and labels
Label = 'label';
Image = 'image';
Fourier = '_Fourier';
FuncSession = '_backward_1';

% .........................................................................Specify Model Fit parameters
ModelFiles = 'Average';

% .........................................................................ROI
ROImg = '3b_hand.nii.gz';
ROI = '3b_hand';

% .........................................................................Switch for functions
% .........................................................................Step One
switch_prepare = false;
switch_copyregister = false;
switch_vol2mat = false;
switch_prepare_resting = false;
switch_resting_vol2mat = false;
switch_BOLD_normalisation = false;
switch_conc_avg = false;
switch_copy = false;
switch_modelfit = false;
switch_mat2vol = false;

%% ........................................................................Preprocessing pipeline
for i_sub=1:size(Subjects, 2)
    
    CurrSubj = Subjects{i_sub};
    CurrSub = Subj{i_sub};
    
    % .....................................................................Make subject directory
    CurrSubjPath = fullfile(RootPath,CurrSubj);
    mkdir (CurrSubjPath);
    
    for i_cond = 1:size(Conditions,2)
        
        CurrCond = Conditions{i_cond};
        
        % .................................................................Make subject functional data directory
        CurrDataPath = fullfile(DataPath,CurrSubj,FldDev,[FldSub CurrSub]);
        
        % .................................................................Make directory and copy files NIFTI
        if switch_prepare
            
            pl_layer_file_preparation(RootPath, DataPath, CurrSubjPath, CurrSubj, CurrSub, CurrCond, Order, RunCon, FldSPM, FldPRF, FldDev, FldSub, ROImg);
            
        end
        
        CurrPathSPM = fullfile(CurrSubjPath,[FldSPM CurrCond]);
        CurrPathPRF = fullfile(CurrSubjPath,[FldPRF CurrCond]);
        
        % .................................................................Copy register.dat file
        
        if switch_copyregister
            
            pl_copyregister(SessPath, CurrPathSPM, CurrSubj, CurrCond, Image, Fourier, FuncSession)
            
        end
        
        % .................................................................Convert Nifti to .mat file
        if switch_vol2mat
            
            pl_vol2mat(CurrSubjPath, FldSPM, CurrCond, Order, RunCon, ROImg, ROI);
            
        end
        
        % .................................................................Copy resting state images
        if switch_prepare_resting
            
            pl_layer_resting_preparation(DataPath, CurrSubjPath, CurrSubj, CurrSub, CurrCond, FldSPM, FldPRF, FldDev, FldSub);
            
        end
        
        % .................................................................Convert Nifti to .mat file
        if switch_resting_vol2mat
            
            pl_resting_vol2mat(CurrSubjPath, FldSPM, CurrCond, ROI);
            
        end
        
        % .................................................................BOLD signal normalisation
        if switch_BOLD_normalisation
            
            pl_BOLD_normalisation(CurrSubjPath, FldSPM, CurrCond, Order, RunCon)
            
        end
        
        % .................................................................Concatenate and average
        
        if switch_conc_avg
            
            cd(CurrPathSPM)
            
            pl_layer_pRF_concatenate(Order, RunCon)
            pl_layer_pRF_average(RunAvg)
            
        end
        
        % .................................................................Copy over ModelFit files
        
        if switch_copy
            
            pl_layer_copyfile(ModelFiles, CurrPathSPM, CurrPathPRF)
            
        end
        
        % .................................................................Gaussian Tuning Curve Vertical
        
        if switch_modelfit
            
            cd(CurrPathPRF)
            
            pl_layer_Gaussian_Tuning_Curve_Vertical(ModelFiles)
            
        end
        
        % .................................................................Convert .mat file to Nifti
        if switch_mat2vol
            
            pl_mat2vol(CurrSubjPath, FldPRF, CurrCond);
            
        end
        
        
    end
    
end

%% ........................................................................Remove SamSrf path
rmpath(genpath('/Users/pliu/Documents/Programmes/samsrf_v9.4'));