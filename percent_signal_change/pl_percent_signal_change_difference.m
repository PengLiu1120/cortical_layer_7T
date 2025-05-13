%% Percent Signal Response Difference
% .........................................................................
% This script compare percent signal change between (D2)+(D3) and (D2+D3).
% .........................................................................
% Inputs
% .........................................................................
% Percent signal change 
% D2
% D3
% D2+D3
% .........................................................................
% Outputs
% .........................................................................
% Difference percent signal response
% .........................................................................
% Written by P.Liu
% Optimized by P.Liu
% Email: peng.liu@med.ovgu.de
% Last modified by P.Liu 10 May 2022
%% ........................................................................Tidy up
clear all
close all
clc

%% ........................................................................Specify parameters
% .........................................................................Specify root path
RootPath = '/Volumes/IKND/AG_Kuehn/Peng/LayerPRF/LayerPRF';

% .........................................................................Specify subjects
Subjects = {'ajz367' 'bkn792' 'bmg520' 'cxc075' 'czg996' 'frj712' 'ggp057' 'gph998' 'gxo876' 'hby152' 'ijt563' 'iwq192' 'kdy341' 'llh150' 'lpr469' 'nhm378' 'oms448' 'qet940' 'qxo538' 'sst050' 'unk742'};
% 'ajz367' 'bkn792' 'bmg520' 'cxc075' 'czg996' 'frj712' 'ggp057' 'gph998' 'gxo876' 'hby152' 'ijt563' 'iwq192' 'kdy341' 'llh150' 'lpr469' 'nhm378' 'oms448' 'qet940' 'qxo538' 'sst050' 'unk742'

% .........................................................................Specify condtions, D2 and D2+D3 or D3 and D2+D3
Conditions = {'D2' 'D3' 'D2+D3'};

% .........................................................................Specify localisers, hand and 3b
Loc = {'3b_hand'};

% .........................................................................Parameters
Results = 'Results';
SignalChange = '%signalchange';

%% ........................................................................Read out percent signal change for each condition
for i_sub=1:size(Subjects, 2)
    
    CurrSubj = Subjects{i_sub};
    
    % .....................................................................Direct to subject path
    CurrSubjPath = fullfile(RootPath,CurrSubj);
    
    for i_loc = 1:size(Loc,2)
        
        CurrLoc = Loc{i_loc};
        
        signal = [];
        
        for i_cond = 1:size(Conditions,2)
            
            CurrCond = Conditions{i_cond};
            
            % .................................................................Direct to % signal change path
            CurrResultPath = fullfile(CurrSubjPath,Results,SignalChange,CurrCond);
            
            cd(CurrResultPath)
            
            signal_change = [CurrLoc '_' CurrCond '_percentsignalchange.mat'];
            
            load(signal_change);
            
            signal = [signal percent_signal_change];
            
        end
        
        differences = signal(:,1) + signal(:,2) - signal(:,3) + 0.5;
        
        % .................................................................Read out coordinates
        
        Label = ['lh.' CurrLoc '_' CurrCond '_percentsignalchange.label'];
        
        vertex = load(Label);
        format long;
        Vertices = vertex(:,1);
        Coordinates_X = vertex(:,2);
        Coordinates_Y = vertex(:,3);
        Coordinates_Z = vertex(:,4);
        
        extra_col = zeros(size(Vertices));
        
        diff_response = [Vertices Coordinates_X Coordinates_Y Coordinates_Z differences extra_col];
        responseName = ['lh.' CurrLoc '_diff_percentsignalchange.label'];
        
        CurrSigPath = fullfile(CurrSubjPath,Results,SignalChange);
        
        cd(CurrSigPath)
            
        dlmwrite(responseName, diff_response, 'delimiter','\t');
        
    end
    
end     