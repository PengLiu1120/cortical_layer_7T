%% Data cleaning
% This script cleans the model filt results based on beta and nR^2
% Data will be cleaned if:
% Beta < 0
% nR^2 too big than 1
%% ........................................................................
% Written by P.Liu
% Optimized by P.Liu
% Email: peng.liu@med.ovgu.de
% Last updated 28 Aug 2022 by P.Liu
%% ........................................................................Tidy up
clear all
close all
clc

%% ........................................................................Get path to SamSrf
addpath(genpath('/Users/pliu/Documents/Programmes/samsrf_v9.4'));

%% ........................................................................Defaults
%..........................................................................Specify modelling path
RootPath = '/Volumes/IKND/AG_Kuehn/Peng/LayerPRF/pRFmodel';

% .........................................................................Specify subjects
Subjects = {'ajz367' 'bkn792' 'bmg520' 'cxc075' 'czg996' 'frj712' 'ggp057' 'gph998' 'gxo876' 'hby152' 'ijt563' 'iwq192' 'kdy341' 'llh150' 'lpr469' 'nhm378' 'oms448' 'qet940' 'qxo538' 'sst050' 'unk742'};

% .........................................................................Specify condtions, D2, D3 and D2+D3
Conditions = {'D2' 'D3' 'D2+D3'};
% 'D2' 'D3' 'D2+D3'

% .........................................................................Parameters
FldPRF = 'pRF_';

%% ........................................................................Data cleaning
for i_sub=1:size(Subjects, 2)
    
    CurrSubj = Subjects{i_sub};
    
    % .....................................................................Subject directory
    CurrSubjPath = fullfile(RootPath,CurrSubj);
    
    for i_cond = 1:size(Conditions,2)
        
        CurrCond = Conditions{i_cond};
        
        % .................................................................PRF directory
        CurrPathPRF = fullfile(CurrSubjPath,[FldPRF CurrCond]);
        
        cd(CurrPathPRF)
        
        % .................................................................Load result .mat file
        load('lh_pTC_ver.mat');
        
        % .................................................................Add new value to Srf.Value
        Values = Srf.Values;
        Values(8,1) = {'orig_nR^2'};
        
        % .................................................................Add original nR^2 to Srf.Data
        Srf.Data(8,:) = Srf.Data(1,:);
        
        % .................................................................Fake R2 vector
        Idx = Srf.Data(4,:) > 0;
        Srf.Data(1,Idx == 0) = 0;
        Srf.Data(1,Idx == 1) = 1;
        
        % .................................................................Data clearning according to Beta
        Srf.Values = Values;
        
        % .................................................................Save cleaned data
        FileName = ['lh_pTC_ver_' CurrCond];
        save(FileName, 'Model', 'Srf', '-v7.3');
        
    end
    
end
%% ........................................................................Remove SamSrf path
rmpath(genpath('/Users/pliu/Documents/Programmes/samsrf_v9.4'));