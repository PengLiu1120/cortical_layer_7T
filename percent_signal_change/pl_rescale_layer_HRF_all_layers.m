%% Rescale layer-specific HRF data
% .........................................................................
% This script reads out the layer-specific HRF data and rescale it between
% 0 and 1 for all individuals, and for young and old group respectively. 
% .........................................................................
% Written by P.Liu
% Optimized by P.Liu
% Email: peng.liu@med.uni-tuebingen.de
% Last modified by P.Liu 13 Jun 2023
%% ........................................................................Tidy up
clear all
close all
clc

%% ........................................................................Set directory
DataDir = '/home/pliu/Documents/LayerPRF/Layer_Specific/08_LayerExtraction/HRF';

%% ........................................................................Set defaults
% .........................................................................Specify subjects
Young = {'frj712' 'gxo876' 'hby152' 'ijt563' 'kdy341' 'lpr469' 'nhm378' 'oms448' 'qet940' 'qxo538' 'unk742'};
Old = {'ajz367' 'bkn792' 'bmg520' 'cxc075' 'czg996' 'ggp057' 'gph998' 'iwq192' 'llh150' 'sst050'};

%% ........................................................................Combine resting state readouts for all young individual
HRF_young = [];

for i_young=1:size(Young, 2)
    
    CurrSubj = Young{i_young};
    
    DataPath = fullfile(DataDir,CurrSubj);
    cd(DataPath);
    
    HRF = [CurrSubj '_HRF_all_layers'];
    load(HRF);
    
    HRF_young = [HRF_young mean_layers];
    
end

cd(DataDir);

Result = 'HRF_young';
save(Result, 'HRF_young');

%% ........................................................................Combine resting state readouts for all old individual
HRF_old = [];

for i_old=1:size(Old, 2)
    
    CurrSubj = Old{i_old};
    
    DataPath = fullfile(DataDir,CurrSubj);
    cd(DataPath);
    
    HRF = [CurrSubj '_HRF_all_layers'];
    load(HRF);
    
    HRF_old = [HRF_old mean_layers];
    
end

cd(DataDir);

Result = 'HRF_old';
save(Result, 'HRF_old');

%% ........................................................................Combine young and old readouts and rescale for both age groups
HRF_all = [HRF_young HRF_old];

Result = 'HRF_all';
save(Result, 'HRF_all');

res_HRF_all = rescale(HRF_all);

Res = 'res_HRF_all';
save(Res, 'res_HRF_all');

res_HRF_young = res_HRF_all(:,1:11);

Res = 'res_HRF_young';
save(Res, 'res_HRF_young');

res_HRF_old = res_HRF_all(:,12:21);

Res = 'res_HRF_old';
save(Res, 'res_HRF_old');