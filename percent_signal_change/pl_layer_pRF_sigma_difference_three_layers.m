%% Layer-specific Percent _sigma Response Difference for Three Layers
% .........................................................................
% This script compare percent _sigma change between (D2)+(D3) and (D2+D3)
% for each layer and each individual.
% .........................................................................
% Inputs
% .........................................................................
% Percent _sigma change 
% D2
% D3
% D2+D3
% .........................................................................
% Outputs
% .........................................................................
% Difference percent _sigma response
% .........................................................................
% Written by P.Liu
% Optimized by P.Liu
% Email: peng.liu@med.uni-tuebingen.de
% Last modified by P.Liu 13 Jun 2023
%% ........................................................................Tidy up
clear all
close all
clc

%% ........................................................................Specify parameters
% .........................................................................Specify root path
DataPath = '/home/pliu/Documents/LayerPRF/Layer_Specific/08_LayerExtraction/pRF';
ResultPath = '/home/pliu/Documents/LayerPRF/Layer_Specific/08_LayerExtraction/pRF/Difference';

% .........................................................................Specify subjects
Young = {'frj712' 'gxo876' 'hby152' 'ijt563' 'kdy341' 'lpr469' 'nhm378' 'oms448' 'qet940' 'qxo538' 'unk742'};
Old = {'ajz367' 'bkn792' 'bmg520' 'cxc075' 'czg996' 'ggp057' 'gph998' 'iwq192' 'llh150' 'sst050'};

% .........................................................................Specify condtions, D2 and D2+D3 or D3 and D2+D3
Conditions = {'D2+D3' 'D2' 'D3'};

%% ........................................................................Read out percent _sigma change for each condition
cd(DataPath);

D2D3_sigma_young = load([Conditions{1} '_young_sigma']);
D2D3_sigma_old = load([Conditions{1} '_old_sigma']);

D2_sigma_young = load([Conditions{2} '_young_sigma']);
D2_sigma_old = load([Conditions{2} '_old_sigma']);

D3_sigma_young = load([Conditions{3} '_young_sigma']);
D3_sigma_old = load([Conditions{3} '_old_sigma']);

sigma_difference_young = [];
sigma_difference_old = [];

for i_young=1:size(Young, 2)
        
        CurrSubj = Young{i_young};
        
        sigma_difference = D2_sigma_young.sigma_all(:,i_young) + D3_sigma_young.sigma_all(:,i_young) - D2D3_sigma_young.sigma_all(:,i_young);
        
        cd(ResultPath);
        Result = [CurrSubj '_sigma_difference'];
        save(Result, 'sigma_difference');
        
        sigma_difference_young = [sigma_difference_young sigma_difference];
        
        cd(DataPath);
        GroupResult = 'sigma_difference_young';
        save(GroupResult, 'sigma_difference_young');
        
end

for i_old=1:size(Old, 2)
        
        CurrSubj = Old{i_old};
        
        sigma_difference = D2_sigma_old.sigma_all(:,i_old) + D3_sigma_old.sigma_all(:,i_old) - D2D3_sigma_old.sigma_all(:,i_old);
        
        cd(ResultPath);
        Result = [CurrSubj '_sigma_difference'];
        save(Result, 'sigma_difference');
        
        sigma_difference_old = [sigma_difference_old sigma_difference];
        
        cd(DataPath);
        GroupResult = 'sigma_difference_old';
        save(GroupResult, 'sigma_difference_old');
        
end