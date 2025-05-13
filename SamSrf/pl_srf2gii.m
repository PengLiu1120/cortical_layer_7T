%% Surface to Gifti conversion
% This script convert SamSrf model fit on surface level to Gifti
% Later then can use mri_surf2vol to convert Gifti to Nifti on volume level
%% ........................................................................
% Written by P.Liu
% Optimized by P.Liu
% Email: peng.liu@med.ovgu.de
% Last updated 15 Sep 2022 by P.Liu
%% ........................................................................Tidy up
clear all
close all
clc

%% ........................................................................Get path to SamSrf
addpath(genpath('/Users/pliu/Documents/Programmes/samsrf_v9.4'));

%% ........................................................................Defaults
%..........................................................................Specify modelling path
RootPath = '/Users/pliu/Documents/LayerPRF/PRFmodelling';

% .........................................................................Specify subjects
Subjects = {'frj712'};
% 'frj712' 'lpr469' 'qet940' 'qxo538' 'oms448' 'hby152' 'ijt563'
% 'ejk164' 'gxo876'

% .........................................................................Specify condtions, D2, D3 and D2+D3
Conditions = {'D2'};
% 'D2' 'D3' 'D2+D3'

% .........................................................................Specify hemisphere, lh and rh
Hemis = {'lh'};
% 'lh' 'rh'

% .........................................................................Parameters
FldPRF = 'pRF_';

%% ........................................................................Convert from SamSrf .mat to Gifti
for i_sub=1:size(Subjects, 2)
    
    CurrSubj = Subjects{i_sub};
    
    % .....................................................................Subject directory
    CurrSubjPath = fullfile(RootPath,CurrSubj);
    
    for i_cond = 1:size(Conditions,2)
        
        CurrCond = Conditions{i_cond};
        
        % .................................................................cd current condition path
        CurrPathPRF = fullfile(CurrSubjPath,[FldPRF CurrCond]);
        cd(CurrPathPRF)
        
        for i_hemi = 1:size(Hemis,2)
            
            CurrHemi = Hemis{i_hemi};
            
            % .............................................................Load result .mat file
            % load('lh_pTC_ver.mat');
            load([CurrHemi '_pTC_ver_' CurrCond]);
            % load([CurrHemi '_pTC_ver']);
            
            % .............................................................Expand Srf
            Srf = samsrf_expand_srf(Srf);
            
            % .............................................................srf2gii for Mu and Sigma
            
            samsrf_srf2gii(Srf, 2)
            samsrf_srf2gii(Srf, 3)
            
        end
        
    end
    
end

%% ........................................................................Remove SamSrf path
rmpath(genpath('/Users/pliu/Documents/Programmes/samsrf_v9.2.2'));