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
Data_Dir = '/home/pliu/Documents/LayerPRF/Layer_Specific/08_LayerExtraction/pRF';

%% ........................................................................Set defaults
% .........................................................................Specify subjects
Subjects = {'ajz367' 'bkn792' 'bmg520' 'cxc075' 'czg996' 'frj712' 'ggp057' 'gph998' 'gxo876' 'hby152' 'ijt563' 'iwq192' 'kdy341' 'llh150' 'lpr469' 'nhm378' 'oms448' 'qet940' 'qxo538' 'sst050' 'unk742'};
Young = {'frj712' 'gxo876' 'hby152' 'ijt563' 'kdy341' 'lpr469' 'nhm378' 'oms448' 'qet940' 'qxo538' 'unk742'};
Old = {'ajz367' 'bkn792' 'bmg520' 'cxc075' 'czg996' 'ggp057' 'gph998' 'iwq192' 'llh150' 'sst050'};

% .........................................................................Specify condtions, D2, D3 and D2+D3
Conditions = {'D2+D3' 'D2' 'D3'};

% .........................................................................Specify layers, Deep, Middle and Superficial
Layers = {'Deep' 'Middle' 'Superficial'};

% .........................................................................Specify age, young and old
Age = {'young' 'old'};

%% ........................................................................Combine pRF sigma results for all individuals
for i_cond = 1:size(Conditions,2)
    
    CurrCond = Conditions{i_cond};
    
    for i_sub=1:size(Subjects, 2)
        
        CurrSubj = Subjects{i_sub};
        
        sigma_individual = [];
        
        Result_Dir = fullfile(Data_Dir,CurrCond,CurrSubj);
        cd(Result_Dir);
        
        load([CurrSubj '_pRF_Sigma_Deep_layers.mat']);
        
        sigma_individual = [sigma_individual; mean_deep(1)];
        
        load([CurrSubj '_pRF_Sigma_Middle_layers.mat']);
        
        sigma_individual = [sigma_individual; mean_middle(1)];
        
        load([CurrSubj '_pRF_Sigma_Superficial_layers.mat']);
        
        sigma_individual = [sigma_individual; mean_superficial(1)];
        
        FileName = [CurrSubj '_' CurrCond '_sigma'];
        save(FileName, 'sigma_individual');
        
    end
    
end

%% ........................................................................Combine pRF sigma results for all young and all old respectively
for i_cond = 1:size(Conditions,2)
    
    CurrCond = Conditions{i_cond};
    
    sigma_all = [];
    
    for i_young=1:size(Young, 2)
        
        CurrSubj = Young{i_young};
        
        Result_Dir = fullfile(Data_Dir,CurrCond,CurrSubj);
        cd(Result_Dir);
        
        load([CurrSubj '_' CurrCond '_sigma.mat']);
        
        sigma_all = [sigma_all sigma_individual];
        
    end
    
    Save_Dir = fullfile(Data_Dir,CurrCond);
    cd(Save_Dir);
    
    FileName = [CurrCond '_young_sigma'];
    save(FileName, 'sigma_all');
    
end
        
for i_cond = 1:size(Conditions,2)
    
    CurrCond = Conditions{i_cond};
    
    sigma_all = [];
    
    for i_old=1:size(Old, 2)
        
        CurrSubj = Old{i_old};
        
        Result_Dir = fullfile(Data_Dir,CurrCond,CurrSubj);
        cd(Result_Dir);
        
        load([CurrSubj '_' CurrCond '_sigma.mat']);
        
        sigma_all = [sigma_all sigma_individual];
        
    end
    
    Save_Dir = fullfile(Data_Dir,CurrCond);
    cd(Save_Dir);
    
    FileName = [CurrCond '_old_sigma'];
    save(FileName, 'sigma_all');
    
end

%% ........................................................................Combine and rescale pRF sigma results for all young and all old separately
sigma = [];

for i_age = 1:size(Age,2)
    
    CurrAge = Age{i_age};
    
    sigma_condition = [];
    
    for i_cond = 1:size(Conditions,2)
    
    CurrCond = Conditions{i_cond};
    
    Result_Dir = fullfile(Data_Dir,CurrCond);
    cd(Result_Dir);
    
    load([CurrCond '_' CurrAge '_sigma.mat']);
    
    sigma_condition = [sigma_condition sigma_all];
    
    end
    
    sigma = [sigma sigma_condition];
    
end

res_sigma = rescale(sigma);

cd(Data_Dir);

res_sigma_young = res_sigma(:,1:11);
Res_all = 'D2+D3_res_sigma_young';
save(Res_all, 'res_sigma_young');

res_sigma_young = res_sigma(:,12:22);
Res_all = 'D2_res_sigma_young';
save(Res_all, 'res_sigma_young');

res_sigma_young = res_sigma(:,23:33);
Res_all = 'D3_res_sigma_young';
save(Res_all, 'res_sigma_young');

res_sigma_old = res_sigma(:,34:43);
Res_all = 'D2+D3_res_sigma_old';
save(Res_all, 'res_sigma_old');

res_sigma_old = res_sigma(:,44:53);
Res_all = 'D2_res_sigma_old';
save(Res_all, 'res_sigma_old');

res_sigma_old = res_sigma(:,54:63);
Res_all = 'D3_res_sigma_old';
save(Res_all, 'res_sigma_old');