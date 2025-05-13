%% Rescale layer-specific percent signal change results
% .........................................................................
% This script reads out the layer-specific percent signal change result and
% rescale it between 0 and 1 for all individuals, and for young and old
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
DataDir = '/home/pliu/Documents/LayerPRF/Layer_Specific/08_LayerExtraction/pRF_Layer';

%% ........................................................................Set defaults
% .........................................................................Specify subjects
Young = {'frj712' 'gxo876' 'hby152' 'ijt563' 'kdy341' 'lpr469' 'nhm378' 'oms448' 'qet940' 'qxo538' 'unk742'};
Old = {'ajz367' 'bkn792' 'bmg520' 'cxc075' 'czg996' 'ggp057' 'gph998' 'iwq192' 'llh150' 'sst050'};

% .........................................................................Specify condtions, D2 and D2+D3 or D3 and D2+D3
Conditions = {'D2+D3' 'D2' 'D3'};

%% ........................................................................Combine percent signal change results for all young individual
sigma_young_all = [];

for i_cond = 1:size(Conditions,2)
    
    CurrCond = Conditions{i_cond};
    
    sigma_young = [];
    
    for i_young=1:size(Young, 2)
        
        CurrSubj = Young{i_young};
        
        DataPath = fullfile(DataDir, CurrCond, CurrSubj);
        cd(DataPath);
        
        sigma = [CurrSubj '_' CurrCond '_all'];
        load(sigma);
        
        sigma_young = [sigma_young mean_layers];
        
    end
    
    cd(DataDir);
    
    Result = [CurrCond '_sigma_young'];
    save(Result, 'sigma_young');
    
    sigma_young_all = [sigma_young_all sigma_young];
    
end

%% ........................................................................Combine percent signal change results for all old individual
sigma_old_all = [];

for i_cond = 1:size(Conditions,2)
    
    CurrCond = Conditions{i_cond};
    
    sigma_old = [];
    
    for i_old=1:size(Old, 2)
        
        CurrSubj = Old{i_old};
        
        DataPath = fullfile(DataDir, CurrCond, CurrSubj);
        cd(DataPath);
        
        sigma = [CurrSubj '_' CurrCond '_all'];
        load(sigma);
        
        sigma_old = [sigma_old mean_layers];
        
    end
    
    cd(DataDir);
    
    Result = [CurrCond '_sigma_old'];
    save(Result, 'sigma_old');
    
    sigma_old_all = [sigma_old_all sigma_old];
    
end

%% ........................................................................Combine young and old results and rescale for both age groups
cd(DataDir);

sigma_all = [sigma_young_all sigma_old_all];

res_sigma_all = rescale(sigma_all);

res_sigma_young = res_sigma_all(:,1:11);
Res_all = 'D2+D3_res_all_sigma_young';
save(Res_all, 'res_sigma_young');

res_sigma_young = res_sigma_all(:,12:22);
Res_all = 'D2_res_all_sigma_young';
save(Res_all, 'res_sigma_young');

res_sigma_young = res_sigma_all(:,23:33);
Res_all = 'D3_res_all_sigma_young';
save(Res_all, 'res_sigma_young');

res_sigma_old = res_sigma_all(:,34:43);
Res_all = 'D2+D3_res_all_sigma_old';
save(Res_all, 'res_sigma_old');

res_sigma_old = res_sigma_all(:,44:53);
Res_all = 'D2_res_all_sigma_old';
save(Res_all, 'res_sigma_old');

res_sigma_old = res_sigma_all(:,54:63);
Res_all = 'D3_res_all_sigma_old';
save(Res_all, 'res_sigma_old');