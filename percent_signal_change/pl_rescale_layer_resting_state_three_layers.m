%% Combine and percent signal change results for three layers
% .........................................................................
% This script combines the percent signal change results for
% three layers (i.e., deep, middle and superficial) for each individual.
% .........................................................................
% Written by P.Liu
% Optimized by P.Liu
% Email: peng.liu@med.uni-tuebingen.de
% Last modified by P.Liu 06 Jun 2023
%% ........................................................................Tidy up
clear all
close all
clc

%% ........................................................................Set directory
Data_Dir = '/home/pliu/Documents/LayerPRF/Layer_Specific/08_LayerExtraction/Resting_State';

%% ........................................................................Set defaults
% .........................................................................Specify subjects
Subjects = {'ajz367' 'bkn792' 'bmg520' 'cxc075' 'czg996' 'frj712' 'ggp057' 'gph998' 'gxo876' 'hby152' 'ijt563' 'iwq192' 'kdy341' 'llh150' 'lpr469' 'nhm378' 'oms448' 'qet940' 'qxo538' 'sst050' 'unk742'};
Young = {'frj712' 'gxo876' 'hby152' 'ijt563' 'kdy341' 'lpr469' 'nhm378' 'oms448' 'qet940' 'qxo538' 'unk742'};
Old = {'ajz367' 'bkn792' 'bmg520' 'cxc075' 'czg996' 'ggp057' 'gph998' 'iwq192' 'llh150' 'sst050'};

% .........................................................................Specify layers, Deep, Middle and Superficial
Layers = {'Deep' 'Middle' 'Superficial'};

% .........................................................................Specify age, young and old
Age = {'young' 'old'};

%% ........................................................................Combine pRF resting_state results for all individuals

for i_sub=1:size(Subjects, 2)
    
    CurrSubj = Subjects{i_sub};
    
    resting_state_individual = [];
    
    Result_Dir = fullfile(Data_Dir,CurrSubj);
    cd(Result_Dir);
    
    load([CurrSubj '_resting_state_Deep_layers_3b.mat']);
    
    resting_state_individual = [resting_state_individual; mean_deep(1)];
    
    load([CurrSubj '_resting_state_Middle_layers_3b.mat']);
    
    resting_state_individual = [resting_state_individual; mean_middle(1)];
    
    load([CurrSubj '_resting_state_Superficial_layers_3b.mat']);
    
    resting_state_individual = [resting_state_individual; mean_superficial(1)];
    
    FileName = [CurrSubj '_resting_state'];
    save(FileName, 'resting_state_individual');
    
end

%% ........................................................................Combine pRF resting_state results for all young and all old respectively
resting_state_young = [];

for i_young=1:size(Young, 2)
    
    CurrSubj = Young{i_young};
    
    Result_Dir = fullfile(Data_Dir,CurrSubj);
    cd(Result_Dir);
    
    load([CurrSubj '_resting_state.mat']);
    
    resting_state_young = [resting_state_young resting_state_individual];
    
end

cd(Data_Dir);

FileName = 'resting_state_young';
save(FileName, 'resting_state_young');

resting_state_old = [];

for i_old=1:size(Old, 2)
    
    CurrSubj = Old{i_old};
    
    Result_Dir = fullfile(Data_Dir,CurrSubj);
    cd(Result_Dir);
    
    load([CurrSubj '_resting_state.mat']);
    
    resting_state_old = [resting_state_old resting_state_individual];
    
end

cd(Data_Dir);

FileName = 'resting_state_old';
save(FileName, 'resting_state_old');

%% ........................................................................Combine and rescale pRF resting_state results for all young and all old separately
resting_state_all = [resting_state_young resting_state_old];

Result = 'resting_state_all';
save(Result, 'resting_state_all');

res_resting_state_all = rescale(resting_state_young);

Res = 'res_resting_state_all';
save(Res, 'res_resting_state_all');

res_resting_state_young = res_resting_state_all(:,1:11);

Res = 'res_resting_state_young';
save(Res, 'res_resting_state_young');

res_resting_state_old = res_resting_state_all(:,12:21);

Res = 'res_resting_state_old';
save(Res, 'res_resting_state_old');