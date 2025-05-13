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
Data_Dir = '/media/pliu/LayerPRF/LayerMapping/08_LayerExtraction/Fourier';

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

%% ........................................................................Combine percent signal change results for all individuals
for i_cond = 1:size(Conditions,2)
    
    CurrCond = Conditions{i_cond};
    
    for i_sub=1:size(Subjects, 2)
        
        CurrSubj = Subjects{i_sub};
        
        percent_signal_change_individual = [];
        
        Result_Dir = fullfile(Data_Dir,CurrCond,CurrSubj);
        cd(Result_Dir);
        
        for i_layer = 1:size(Layers,2)
            
            CurrLayer = Layers{i_layer};
            
            load([CurrCond '_' CurrLayer '_percentsignalchange.mat']);
            
            percent_signal_change_individual = [percent_signal_change_individual; percent_signal_change];
            
        end
        
        FileName = [CurrSubj '_' CurrCond '_percent_signal_change'];
        save(FileName, 'percent_signal_change_individual');
        
    end
    
end

%% ........................................................................Combine percent signal change results for all young and all old respectively
for i_cond = 1:size(Conditions,2)
    
    CurrCond = Conditions{i_cond};
    
    percent_signal_change_all = [];
    
    for i_young=1:size(Young, 2)
        
        CurrSubj = Young{i_young};
        
        Result_Dir = fullfile(Data_Dir,CurrCond,CurrSubj);
        cd(Result_Dir);
        
        load([CurrSubj '_' CurrCond '_percent_signal_change.mat']);
        
        percent_signal_change_all = [percent_signal_change_all percent_signal_change_individual];
        
    end
    
    Save_Dir = fullfile(Data_Dir,CurrCond);
    cd(Save_Dir);
    
    FileName = [CurrCond '_young_percent_signal_change'];
    save(FileName, 'percent_signal_change_all');
    
end
        
for i_cond = 1:size(Conditions,2)
    
    CurrCond = Conditions{i_cond};
    
    percent_signal_change_all = [];
    
    for i_old=1:size(Old, 2)
        
        CurrSubj = Old{i_old};
        
        Result_Dir = fullfile(Data_Dir,CurrCond,CurrSubj);
        cd(Result_Dir);
        
        load([CurrSubj '_' CurrCond '_percent_signal_change.mat']);
        
        percent_signal_change_all = [percent_signal_change_all percent_signal_change_individual];
        
    end
    
    Save_Dir = fullfile(Data_Dir,CurrCond);
    cd(Save_Dir);
    
    FileName = [CurrCond '_old_percent_signal_change'];
    save(FileName, 'percent_signal_change_all');
    
end

%% ........................................................................Combine and rescale percent signal change results for all young and all old separately
percent_signal_change = [];

for i_age = 1:size(Age,2)
    
    CurrAge = Age{i_age};
    
    percent_signal_change_condition = [];
    
    for i_cond = 1:size(Conditions,2)
    
    CurrCond = Conditions{i_cond};
    
    Result_Dir = fullfile(Data_Dir,CurrCond);
    cd(Result_Dir);
    
    load([CurrCond '_' CurrAge '_percent_signal_change.mat']);
    
    percent_signal_change_condition = [percent_signal_change_condition percent_signal_change_all];
    
    end
    
    percent_signal_change = [percent_signal_change percent_signal_change_condition];
    
end

res_percent_signal_change = rescale(percent_signal_change);

cd(Data_Dir);

res_percent_signal_change_young = res_percent_signal_change(:,1:11);
Res_all = 'D2+D3_res_percent_signal_change_young';
save(Res_all, 'res_percent_signal_change_young');

res_percent_signal_change_young = res_percent_signal_change(:,12:22);
Res_all = 'D2_res_percent_signal_change_young';
save(Res_all, 'res_percent_signal_change_young');

res_percent_signal_change_young = res_percent_signal_change(:,23:33);
Res_all = 'D3_res_percent_signal_change_young';
save(Res_all, 'res_percent_signal_change_young');

res_percent_signal_change_old = res_percent_signal_change(:,34:43);
Res_all = 'D2+D3_res_percent_signal_change_old';
save(Res_all, 'res_percent_signal_change_old');

res_percent_signal_change_old = res_percent_signal_change(:,44:53);
Res_all = 'D2_res_percent_signal_change_old';
save(Res_all, 'res_percent_signal_change_old');

res_percent_signal_change_old = res_percent_signal_change(:,54:63);
Res_all = 'D3_res_percent_signal_change_old';
save(Res_all, 'res_percent_signal_change_old');