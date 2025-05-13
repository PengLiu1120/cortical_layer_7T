%% Plot timecourses for vertices
% This script plot timecourses for selected vertices, compare fitted
% timecourses with observed ones using samsrf_fitvsobs
% .........................................................................

% .........................................................................
% Written by P.Liu
% Optimized by P.Liu
% Email: peng.liu@med.ovgu.de
% Last updated 20 May 2022 by P.Liu
%% ........................................................................Tidy up
clear all
close all
clc

%% ........................................................................Get path to SamSrf
addpath(genpath('/Users/pliu/Documents/Programmes/samsrf_v8.0'));

%% ........................................................................Defaults
%..........................................................................Specify modelling path
RootPath = '/Users/pliu/Documents/LayerPRF/PRFmodelling';

% .........................................................................Specify subjects
Subjects = {'frj712'};
% 'frj712' 'lpr469' 'qet940' 'qxo538' 'oms448'

% .........................................................................Specify condtions, D2, D3 and D2+D3
Conditions = {'D2+D3'};
% 'D2' 'D3' 'D2+D3'

% .........................................................................Specify vertices
v = {'54958', '59789', '37499', '37501'};

% .........................................................................Parameters
FldPRF = 'pRF_';
FldFig = 'figures';

%% ........................................................................Plot timecourse for seleted vertex in finefit
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
        cd(CurrPathPRF)
        
        % .................................................................Load result .mat file
        load('lh_pTC_ver.mat');
        
        % .................................................................Expand Srf
        Srf = samsrf_expand_srf(Srf);
        
        for i_v = 1:size(v, 2)
            
            Curr_v = v{i_v};
            
            vertex = str2double(Curr_v);
        
            [S, X, Y] = samsrf_fitvsobs(Srf, Model, vertex);
            
            cd(CurrFigPath)
            
            FigName = ['FineFit_Timecourse_' Curr_v '_' CurrCond '.png'];
            
            saveas(gcf, FigName)
            
            close all
            
        end
        
    end
    
end
%% ........................................................................Plot timecourse for seleted vertex in coarse fit
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
        cd(CurrPathPRF)
        
        % .................................................................Load result .mat file
        load('lh_pTC_ver_CrsFit.mat');
        
        % .................................................................Expand Srf
        Srf = samsrf_expand_srf(Srf);
        
        for i_v = 1:size(v, 2)
            
            Curr_v = v{i_v};
            
            vertex = str2double(Curr_v);
        
            [S, X, Y] = samsrf_fitvsobs(Srf, Model, vertex);
            
            cd(CurrFigPath)
            
            FigName = ['CoarseFit_Timecourse_' Curr_v '_' CurrCond '.png'];
            
            saveas(gcf, FigName)
            
            close all
            
        end
        
    end
    
end