%% Extract layer-specific pRF size (sigma) results
% .........................................................................
% This script extracts the rescaled layer-specific pRF size (sigma) results
% for young and old group respectively, and for each individual.
% .........................................................................
% Written by P.Liu
% Optimized by P.Liu
% Email: peng.liu@med.ovgu.de
% Last modified by P.Liu 30 May 2023
%% ........................................................................Tidy up
clear all
close all
clc

%% ........................................................................Set directory
Result_Dir = '/media/pliu/LayerPRF/LayerMapping/08_LayerExtraction/pRF_Layer';
cd(Result_Dir);

FigureAllPath = '/home/pliu/Documents/LayerPRF/Figures';
FigurePath = '/home/pliu/Documents/LayerPRF/Figures/pRF_size_layer_individual';

%% ........................................................................Set defaults
% .........................................................................Specify subjects
Young = {'frj712' 'gxo876' 'hby152' 'ijt563' 'kdy341' 'lpr469' 'nhm378' 'oms448' 'qet940' 'qxo538' 'unk742'};
Old = {'ajz367' 'bkn792' 'bmg520' 'cxc075' 'czg996' 'ggp057' 'gph998' 'iwq192' 'llh150' 'sst050'};

% .........................................................................Specify condtions, D2 and D2+D3 or D3 and D2+D3
Conditions = {'D2+D3' 'D2' 'D3'};

%% ........................................................................Combine pRF size (sigma) results for all young individual
for i_cond = 1:size(Conditions,2)
    
    CurrCond = Conditions{i_cond};
    
    sigma_young = [];
    
    for i_young=1:size(Young, 2)
        
        CurrSubj = Young{i_young};
        
        DataPath = fullfile(Result_Dir,CurrCond,CurrSubj);
        cd(DataPath);
        
        Sigma = [CurrSubj '_' CurrCond '_all'];
        load(Sigma);
        
        sigma_young = [sigma_young mean_layers];
        
    end
    
    cd(Result_Dir);
    
    Result = [CurrCond '_sigma_young'];
    save(Result, 'sigma_young');
    
end

for i_cond = 1:size(Conditions,2)
    
    CurrCond = Conditions{i_cond};
    
    sigma_old = [];
    
    for i_old=1:size(Old, 2)
        
        CurrSubj = Old{i_old};
        
        DataPath = fullfile(Result_Dir,CurrCond,CurrSubj);
        cd(DataPath);
        
        Sigma = [CurrSubj '_' CurrCond '_all'];
        load(Sigma);
        
        sigma_old = [sigma_old mean_layers];
        
    end
    
    cd(Result_Dir);
    
    Result = [CurrCond '_sigma_old'];
    save(Result, 'sigma_old');
    
end