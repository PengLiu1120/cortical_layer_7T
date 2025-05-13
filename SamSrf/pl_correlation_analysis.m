%% Correlation analysis
% This script perform correlation analysis using function samsrf_cfcorr
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

%% ........................................................................Defaults
%..........................................................................Specify modelling path
RootPath = '/Volumes/IKND/AG_Kuehn/Peng/LayerPRF/pRFmodel';

% .........................................................................Specify subjects
Subjects = {'ajz367' 'bkn792' 'bmg520' 'cxc075' 'czg996' 'ejk164' 'frj712' 'ggp057' 'gph998' 'gxo876' 'hby152' 'ijt563' 'iwq192' 'kdy341' 'llh150' 'lpr469' 'nhm378' 'oms448' 'qet940' 'qxo538' 'sst050' 'unk742'};
% .........................................................................Specify condtions, D2, D3 and D2+D3
Conditions = {'D2' 'D3' 'D2+D3'};

% .........................................................................Specify ROIs, 3b hand and BA3b
ROI = {'3b_hand'};

% .........................................................................Parameters
FldPRF = 'pRF_';
FldFig = 'figures';

%% ........................................................................Correlation analysis
for i_sub=1:size(Subjects, 2)
    
    CurrSubj = Subjects{i_sub};
    
    % .....................................................................Subject directory
    CurrSubjPath = fullfile(RootPath,CurrSubj);
    
    % .....................................................................Figure directory
    CurrFigPath = fullfile(RootPath,CurrSubj,FldFig);
    mkdir (CurrFigPath);
    
    for i_cond = 1:size(Conditions,2)
        
        CurrCond = Conditions{i_cond};
        
        for i_ROI = 1:size(ROI,2)
            
            CurrROI = ROI{i_ROI};
            
            Fig = sprintf([CurrROI '_Crs_Figure%d.png'], i_cond);
            
            % .............................................................PRF directory
            CurrPathPRF = fullfile(CurrSubjPath,[FldPRF CurrCond]);
            
            cd(CurrPathPRF)
            
            % .............................................................Load S and X (predicted timecourse)
            load('src_pTC_ver.mat');
            
            % .............................................................Load model and Y (observed timecourse)
            load('lh_pTC_ver_CrsFit.mat');
            
            Srf = samsrf_expand_srf(Srf);
            
            Loc = samsrf_loadlabel(['lh.' CurrROI '_Loc']);
            
            Y = Srf.Y(:,Loc);
            
            cd(CurrFigPath)
            
            everynthplots = 1:10:size(Y,2);
            nplot = 1;
            ncolumn = ceil(size(everynthplots,2)/4);
            nrow = ceil(size(everynthplots,2)/ncolumn);
            
            figure('Units', 'Normalized', 'Position', [0 0 1 1]);
            
            for i_plot = everynthplots
                
                subplot(ncolumn,nrow,nplot)
                
                R = pl_cfcorr(Y(:,i_plot), X, S, Model);
                
                nplot = nplot + 1;
                
                set(gca,'fontsize', 18);
                
            end
            
            saveas(gcf, Fig);
            
            close all
            
        end
        
    end
    
end
%% ........................................................................Remove SamSrf path
rmpath(genpath('/Users/pliu/Documents/Programmes/samsrf_v9.4'));