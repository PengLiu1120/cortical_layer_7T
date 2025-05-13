%% Rescale layer-specific resting state functional data
% .........................................................................
% This script reads out the layer-specific resting state functional data
% and rescale it between 0 and 1 for all individuals, and for young and old
% group respectively.
% .........................................................................
% Written by P.Liu
% Optimized by P.Liu
% Email: peng.liu@med.ovgu.de
% Last modified by P.Liu 13 May 2023
%% ........................................................................Tidy up
clear all
close all
clc

%% ........................................................................Set directory
DataDir = 'D:\LayerPRF\Layer_Specific\08_LayerExtraction\Resting_State';

%% ........................................................................Set defaults
% .........................................................................Specify subjects
Young = {'frj712' 'gxo876' 'hby152' 'ijt563' 'kdy341' 'lpr469' 'nhm378' 'oms448' 'qet940' 'qxo538' 'unk742'};
Old = {'ajz367' 'bkn792' 'bmg520' 'cxc075' 'czg996' 'ggp057' 'gph998' 'iwq192' 'llh150' 'sst050'};

%% ........................................................................Combine resting state readouts for all young individual
resting_state_young = [];

for i_young=1:size(Young, 2)
    
    CurrSubj = Young{i_young};
    
    DataPath = fullfile(DataDir,CurrSubj);
    cd(DataPath);
    
    resting_state = [CurrSubj '_resting_state_all_layers_3b'];
    load(resting_state);
    
    resting_state_young = [resting_state_young mean_layers];
    
end

cd(DataDir);

Result = 'resting_state_young_orig';
save(Result, 'resting_state_young');

%% ........................................................................Combine resting state readouts for all old individual
resting_state_old = [];

for i_old=1:size(Old, 2)
    
    CurrSubj = Old{i_old};
    
    DataPath = fullfile(RootDir, DataDir,CurrSubj);
    cd(DataPath);
    
    resting_state = [CurrSubj '_resting_state_all_layers_3b'];
    load(resting_state);
    
    resting_state_old = [resting_state_old mean_layers];
    
end

cd(DataDir);

Result = 'resting_state_old_orig';
save(Result, 'resting_state_old');

%% ........................................................................Combine young and old readouts and rescale for both age groups
resting_state_all = [resting_state_young resting_state_old];

Result = 'resting_state_all';
save(Result, 'resting_state_all');

res_resting_state_all = rescale(resting_state_all);

Res = 'res_resting_state_all';
save(Res, 'res_resting_state_all');

res_resting_state_young = res_resting_state_all(:,1:11);

Res = 'res_resting_state_young';
save(Res, 'res_resting_state_young');

res_resting_state_old = res_resting_state_all(:,12:21);

Res = 'res_resting_state_old';
save(Res, 'res_resting_state_old');