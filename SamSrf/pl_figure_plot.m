%% Figure plot for SamSrf pRF mapping
% This script plots figures automatically using DiplayMaps tool
% function samsrf_surf
% samsrf_surf(Srf, Mesh, Thrsh, Paths, CamView, MapType, PatchHandle)
% .........................................................................
% Inputs
% .........................................................................
% Srf: Load from *_pTC_ver.mat, need to expand using samsrf_expand_srf
% Mesh: 'inflated'
% Threshold: [0, 0, 1.5, 0, Inf, 0]
% Paths: 3b Hand label and BA3b label
% CamView: [-94 15 1.75]
% MapType: Mu (vertical) and Sigma
% PatchHandle: []
%..........................................................................
% Outputs
% .........................................................................
% .png
% .........................................................................
% Written by P.Liu
% Optimized by P.Liu
% Email: peng.liu@med.ovgu.de
% Last Updated 11 Apr 2022 by P.Liu
%% ........................................................................Tidy up
clear all
close all
clc

%% ........................................................................Get path to SamSrf
addpath(genpath('/Users/pliu/Documents/Programmes/samsrf_v9.4'));

%% ........................................................................Defaults
% .........................................................................Specify Paths
RootPath = '/Volumes/IKND/AG_Kuehn/Peng/LayerPRF/pRFmodel';

% .........................................................................Specify subjects
Subjects = {'ajz367' 'bkn792' 'bmg520' 'cxc075' 'czg996' 'ejk164' 'frj712' 'ggp057' 'gph998' 'gxo876' 'hby152' 'ijt563' 'iwq192' 'kdy341' 'llh150' 'lpr469' 'nhm378' 'oms448' 'qet940' 'qxo538' 'sst050' 'unk742'};

% .........................................................................Specify condtions, D2, D3 and D2+D3
Conditions = {'D2' 'D3' 'D2+D3'};
% 'D2' 'D3' 'D2+D3'

% .........................................................................Specify localisers, hand and 3b
Loc = {'3b_hand'};

%% ........................................................................Specify plot parameters
% .........................................................................Mesh
Mesh = 'inflated';

% .........................................................................Thresholds
Threshold = [0, 0, 1.25, 0, Inf, 0];
% Threshold = [0, 0, 0.6, 0, Inf, 0];

% .........................................................................CamView
CamView =[-94 15 1.75];

% .........................................................................MapType
MapType = {'Mu_v' 'Sigma'};
% 'Mu_v' 'Sigma'
% 'Noise Ceiling'

% .........................................................................Fields
FldSPM = 'spm_';
FldPRF = 'pRF_';
% FldFig = 'figures';
FldFig = 'figures_clean';

%% ........................................................................Plotting all vertices and draw labels
for i_sub=1:size(Subjects, 2)
    
    CurrSubj = Subjects{i_sub};
    
    % .....................................................................Subject directory
    CurrSubjPath = fullfile(RootPath,CurrSubj);
    
    % .....................................................................Figure directory
    CurrFigPath = fullfile(RootPath,CurrSubj,FldFig);
    % mkdir (CurrFigPath);
    
    for i_cond = 1:size(Conditions,2)
        
        CurrCond = Conditions{i_cond};
        
        % .................................................................cd current condition path
        CurrPathPRF = fullfile(CurrSubjPath,[FldPRF CurrCond]);
        cd(CurrPathPRF)
        
        % .................................................................Load result .mat file
        % load('lh_pTC_ver.mat');
        load(['lh_pTC_ver_' CurrCond]);
        
        % .................................................................Expand Srf
        Srf = samsrf_expand_srf(Srf);
        
        for i_loc = 1:size(Loc,2)
            
            CurrLoc = Loc{i_loc};
            
            % .............................................................Specify label
            Label = ['lh.' CurrLoc '_Loc.label'];
            
            % .............................................................Paths
            Paths = {Label; [NaN 1 1 1]};
            
            for i_map = 1:size(MapType,2)
                
                CurrMap = MapType{i_map};
                
                cd(CurrPathPRF);
                
                % .........................................................Surf figure plot
                samsrf_surf(Srf, Mesh, Threshold, Paths, CamView, CurrMap)
                
                cd(CurrFigPath)
                
                % FigName = [CurrMap '_' CurrLoc '_' CurrCond '.png'];
                FigName = [CurrMap '_' CurrLoc '_' CurrCond '_clean.png'];
                
                saveas(gcf, FigName)
                
                close all
                
            end
            
        end
        
    end
    
end

%% ........................................................................Plot vertices in the labels only and draw labels
for i_sub=1:size(Subjects, 2)
    
    CurrSubj = Subjects{i_sub};
    
    % .....................................................................Subject directory
    CurrSubjPath = fullfile(RootPath,CurrSubj);
    
    % .....................................................................Figure directory
    CurrFigPath = fullfile(RootPath,CurrSubj,FldFig);
    
    for i_cond = 1:size(Conditions,2)
        
        CurrCond = Conditions{i_cond};
        
        % .................................................................cd current condition path
        CurrPathPRF = fullfile(CurrSubjPath,[FldPRF CurrCond]);
        
        for i_loc = 1:size(Loc,2)
            
            CurrLoc = Loc{i_loc};
            
            cd(CurrPathPRF)
            
            % .............................................................Load result .mat file
            % load('lh_pTC_ver.mat');
            load(['lh_pTC_ver_' CurrCond]);
            
            % .............................................................Expand Srf
            Srf = samsrf_expand_srf(Srf);
            
            % .............................................................Specify label
            Label = ['lh.' CurrLoc '_Loc.label'];
            Local = samsrf_loadlabel(['lh.' CurrLoc '_Loc']);
            
            for index = 1:size(Srf.Data,2)
                
                if index ~= Local
                    
                    Srf.Data(:,index) = NaN;
                    
                end
                
            end
            
            % .............................................................Paths
            Paths = {Label; [NaN 1 1 1]};
            
            for i_map = 1:size(MapType,2)
                
                CurrMap = MapType{i_map};
                
                cd(CurrPathPRF);
                
                % .........................................................Surf figure plot
                samsrf_surf(Srf, Mesh, Threshold, '', CamView, CurrMap)
                
                cd(CurrFigPath)
                
                % FigName = [CurrMap '_' CurrLoc '_only_' CurrCond '.png'];
                FigName = [CurrMap '_' CurrLoc '_only_' CurrCond '_clean.png'];
                
                saveas(gcf, FigName)
                
                close all
                
            end
            
        end
        
    end
    
end

%% ........................................................................Remove SamSrf path
rmpath(genpath('/Users/pliu/Documents/Programmes/samsrf_v9.4'));