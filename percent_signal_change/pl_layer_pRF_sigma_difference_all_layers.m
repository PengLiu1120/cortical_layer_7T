%% Layer-specific pRF size (sigma) Difference for Each Layer
% .........................................................................
% This script compare pRF size (sigma) between (D2)+(D3) and (D2+D3) for
% each layer and each individual. 
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
DataPath = '/media/pliu/LayerPRF/LayerMapping/08_LayerExtraction/pRF_Layer';
ResultPath = '/media/pliu/LayerPRF/LayerMapping/08_LayerExtraction/pRF_Layer/Difference';

% .........................................................................Specify subjects
Young = {'frj712' 'gxo876' 'hby152' 'ijt563' 'kdy341' 'lpr469' 'nhm378' 'oms448' 'qet940' 'qxo538' 'unk742'};
Old = {'ajz367' 'bkn792' 'bmg520' 'cxc075' 'czg996' 'ggp057' 'gph998' 'iwq192' 'llh150' 'sst050'};

% .........................................................................Specify condtions, D2 and D2+D3 or D3 and D2+D3
Conditions = {'D2+D3' 'D2' 'D3'};

%% ........................................................................Read out percent signal change for each condition
cd(DataPath);

D2D3Sigma_young = load([Conditions{1} '_sigma_young']);
D2D3Sigma_old = load([Conditions{1} '_sigma_old']);

D2Sigma_young = load([Conditions{2} '_sigma_young']);
D2Sigma_old = load([Conditions{2} '_sigma_old']);

D3Sigma_young = load([Conditions{3} '_sigma_young']);
D3Sigma_old = load([Conditions{3} '_sigma_old']);

sigma_difference_young = [];
sigma_difference_old = [];

for i_young=1:size(Young, 2)
        
        CurrSubj = Young{i_young};
        
        sigma_difference = D2Sigma_young.sigma_young(:,i_young) + D3Sigma_young.sigma_young(:,i_young) - D2D3Sigma_young.sigma_young(:,i_young);
        
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
        
        sigma_difference = D2Sigma_old.sigma_old(:,i_old) + D3Sigma_old.sigma_old(:,i_old) - D2D3Sigma_old.sigma_old(:,i_old);
        
        cd(ResultPath);
        Result = [CurrSubj '_sigma_difference'];
        save(Result, 'sigma_difference');
        
        sigma_difference_old = [sigma_difference_old sigma_difference];
        
        cd(DataPath);
        GroupResult = 'sigma_difference_old';
        save(GroupResult, 'sigma_difference_old');
        
end