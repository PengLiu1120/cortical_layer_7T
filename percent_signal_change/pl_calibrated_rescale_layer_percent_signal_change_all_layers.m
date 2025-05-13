%% Calibrate rescaled layer-specific percent signal change results
% .........................................................................
% This script calibrates the rescaled layer-specific percent signal change
% results for each individual
% .........................................................................
% Written by P.Liu
% Optimized by P.Liu
% Email: peng.liu@med.ovgu.de
% Last modified by P.Liu 13 May 2023
%% ........................................................................Tidy up
clear all
close all
clc

%% ........................................................................General Specifications
% .........................................................................Specify RootDir
RootDir = '/media/pliu/LayerPRF/LayerMapping/08_LayerExtraction/Fourier';

% .........................................................................Layer file folder
DataDir = '08_LayerExtraction/Fourier_Layer';
RestingDir = '08_LayerExtraction/Resting_State_mean';

%% ........................................................................Set defaults
% .........................................................................Specify subjects
Young = {'frj712' 'gxo876' 'hby152' 'ijt563' 'kdy341' 'lpr469' 'nhm378' 'oms448' 'qet940' 'qxo538' 'unk742'};
Old = {'ajz367' 'bkn792' 'bmg520' 'cxc075' 'czg996' 'ggp057' 'gph998' 'iwq192' 'llh150' 'sst050'};

% .........................................................................Specify condtions, D2 and D2+D3 or D3 and D2+D3
Conditions = {'D2+D3' 'D2' 'D3'};

%% ........................................................................Calibrate rescaled percent signal change results
% .........................................................................Specify resting state readouts
RestingPath = fullfile(RootDir, RestingDir);
cd(RestingPath);

resting_state_young = 'res_resting_state_young.mat';
load(resting_state_young);
resting_state_old = 'res_resting_state_old.mat';
load(resting_state_old);

% .........................................................................Calculate calibration
percent_signal_change_young_all_calc = [];
percent_signal_change_old_all_calc = [];

for i_cond = 1:size(Conditions,2)
    
    CurrCond = Conditions{i_cond};
    
    ResultPath = fullfile(RootDir, DataDir);
    cd(ResultPath);
    
    percent_signal_change_young = [CurrCond '_res_all_percent_signal_change_young'];
    load(percent_signal_change_young);
    percent_signal_change_old = [CurrCond '_res_all_percent_signal_change_old'];
    load(percent_signal_change_old);
    
    percent_signal_change_young_calc = res_percent_signal_change_young(:,:)-res_resting_state_young(:,:);
    percent_signal_change_old_calc = res_percent_signal_change_old(:,:)-res_resting_state_old(:,:);
    
    Calc_young = [CurrCond '_percent_signal_change_young_calc'];
    save(Calc_young, 'percent_signal_change_young_calc');
    
    Calc_old = [CurrCond '_percent_signal_change_old_calc'];
    save(Calc_old, 'percent_signal_change_old_calc');
    
    percent_signal_change_young_all_calc = [percent_signal_change_young_all_calc percent_signal_change_young_calc];
    percent_signal_change_old_all_calc = [percent_signal_change_old_all_calc percent_signal_change_old_calc];
    
end

%% ........................................................................Rescale the calibrated percent signal change results
percent_signal_change_all_calc = [percent_signal_change_young_all_calc percent_signal_change_old_all_calc];

res_percent_signal_change_all_calc = rescale(percent_signal_change_all_calc);
    
res_percent_signal_change_young = res_percent_signal_change_all_calc(:,1:11);
Res_all = 'D2+D3_calibrated_res_all_percent_signal_change_young';
save(Res_all, 'res_percent_signal_change_young');

res_percent_signal_change_young = res_percent_signal_change_all_calc(:,12:22);
Res_all = 'D2_calibrated_res_all_percent_signal_change_young';
save(Res_all, 'res_percent_signal_change_young');

res_percent_signal_change_young = res_percent_signal_change_all_calc(:,23:33);
Res_all = 'D3_calibrated_res_all_percent_signal_change_young';
save(Res_all, 'res_percent_signal_change_young');

res_percent_signal_change_old = res_percent_signal_change_all_calc(:,34:43);
Res_all = 'D2+D3_calibrated_res_all_percent_signal_change_old';
save(Res_all, 'res_percent_signal_change_old');

res_percent_signal_change_old = res_percent_signal_change_all_calc(:,44:53);
Res_all = 'D2_calibrated_res_all_percent_signal_change_old';
save(Res_all, 'res_percent_signal_change_old');

res_percent_signal_change_old = res_percent_signal_change_all_calc(:,54:63);
Res_all = 'D3_calibrated_res_all_percent_signal_change_old';
save(Res_all, 'res_percent_signal_change_old');