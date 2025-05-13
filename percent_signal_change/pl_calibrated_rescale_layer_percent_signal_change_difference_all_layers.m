%% Layer-specific Percent Signal Response Difference for Each Layer
% .........................................................................
% This script compare percent signal change between (D2)+(D3) and (D2+D3)
% for each layer and each individual.
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
% Last modified by P.Liu 14 May 2023
%% ........................................................................Tidy up
clear all
close all
clc

%% ........................................................................Specify parameters
% .........................................................................Specify root path
DataPath = '/home/pliu/Documents/LayerPRF/LayerMapping/08_LayerExtraction/Fourier_Layer';
ResultPath = '/home/pliu/Documents/LayerPRF/LayerMapping/08_LayerExtraction/Fourier_Layer/Difference';

% .........................................................................Specify subjects
Young = {'frj712' 'gxo876' 'hby152' 'ijt563' 'kdy341' 'lpr469' 'nhm378' 'oms448' 'qet940' 'qxo538' 'unk742'};
Old = {'ajz367' 'bkn792' 'bmg520' 'cxc075' 'czg996' 'ggp057' 'gph998' 'iwq192' 'llh150' 'sst050'};

% .........................................................................Specify condtions, D2 and D2+D3 or D3 and D2+D3
Conditions = {'D2+D3' 'D2' 'D3'};

%% ........................................................................Read out percent signal change for each condition
cd(DataPath);

D2D3Signal_young = load([Conditions{1} '_calibrated_res_all_percent_signal_change_young']);
D2D3Signal_old = load([Conditions{1} '_calibrated_res_all_percent_signal_change_old']);

D2Signal_young = load([Conditions{2} '_calibrated_res_all_percent_signal_change_young']);
D2Signal_old = load([Conditions{2} '_calibrated_res_all_percent_signal_change_old']);

D3Signal_young = load([Conditions{3} '_calibrated_res_all_percent_signal_change_young']);
D3Signal_old = load([Conditions{3} '_calibrated_res_all_percent_signal_change_old']);

res_percent_signal_change_difference_young = [];
res_percent_signal_change_difference_old = [];

for i_young=1:size(Young, 2)
        
        CurrSubj = Young{i_young};
        
        res_percent_signal_change_difference = D2Signal_young.res_percent_signal_change_young(:,i_young) + D3Signal_young.res_percent_signal_change_young(:,i_young) - D2D3Signal_young.res_percent_signal_change_young(:,i_young);
        
        cd(ResultPath);
        Result = [CurrSubj '_calibrated_res_percent_signal_change_difference'];
        save(Result, 'res_percent_signal_change_difference');
        
        res_percent_signal_change_difference_young = [res_percent_signal_change_difference_young res_percent_signal_change_difference];
        
        cd(DataPath);
        GroupResult = 'calibrated_res_all_percent_signal_change_difference_young';
        save(GroupResult, 'res_percent_signal_change_difference_young');
        
end

for i_old=1:size(Old, 2)
        
        CurrSubj = Old{i_old};
        
        res_percent_signal_change_difference = D2Signal_old.res_percent_signal_change_old(:,i_old) + D3Signal_old.res_percent_signal_change_old(:,i_old) - D2D3Signal_old.res_percent_signal_change_old(:,i_old);
        
        cd(ResultPath);
        Result = [CurrSubj '_calibrated_res_percent_signal_change_difference'];
        save(Result, 'res_percent_signal_change_difference');
        
        res_percent_signal_change_difference_old = [res_percent_signal_change_difference_old res_percent_signal_change_difference];
        
        cd(DataPath);
        GroupResult = 'calibrated_res_all_percent_signal_change_difference_old';
        save(GroupResult, 'res_percent_signal_change_difference_old');
        
end