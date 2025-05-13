%% Plot rescaled layer-specific percent signal change results
% .........................................................................
% This script plots the rescaled layer-specific percent signal change
% results for young and old group respectively, and for each individual.
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
Result_Dir = '/home/pliu/Documents/LayerPRF/Layer_Specific/08_LayerExtraction/pRF_Layer';
cd(Result_Dir);

FigurePath = '/home/pliu/Documents/LayerPRF/Figures/pRF_size_layers';

%% ........................................................................Set defaults
% .........................................................................Specify subjects
Young = {'frj712' 'gxo876' 'hby152' 'ijt563' 'kdy341' 'lpr469' 'nhm378' 'oms448' 'qet940' 'qxo538' 'unk742'};
Old = {'ajz367' 'bkn792' 'bmg520' 'cxc075' 'czg996' 'ggp057' 'gph998' 'iwq192' 'llh150' 'sst050'};

% .........................................................................Specify condtions, D2, D3 and D2+D3
Conditions = {'D2+D3' 'D2' 'D3'};

% .........................................................................Specify rescaled layer-specific percent signal change results
sigma_young='D2+D3_res_all_sigma_young.mat';
load(sigma_young);
sigma_old='D2+D3_res_all_sigma_old.mat';
load(sigma_old);

% .........................................................................Specify plot colours
colorspec_group = {'0000e1','e10000','868584'};

color_young = sscanf(colorspec_group{1},'%2x%2x%2x',[1 3])/255;
color_old = sscanf(colorspec_group{2},'%2x%2x%2x',[1 3])/255;
color_ind = sscanf(colorspec_group{3},'%2x%2x%2x',[1 3])/255;

% .........................................................................Specify group numbers
n_young = 11;
n_old = 10;

%% ........................................................................Calculate mean
mean_res_young = nanmean(res_sigma_young,2);
mean_res_old = nanmean(res_sigma_old,2);

%% ........................................................................Plot figures for young
figure(1)

set(gca,'xtick',[],'ytick',[],'LineWidth',4,'fontsize',12,'FontWeight','bold');
set(gca,'Layer', 'top');

ylim([0 1])
hold on;

for ind = 1:n_young    
    figure(1)        
    plot(res_sigma_young(:,ind),'Color',color_ind,'LineWidth',4)
    hold on;
end
 
figure(1);
plot(mean_res_young,'Color',color_young,'LineWidth',6);

%% ........................................................................Plot figures for old
figure(2)

set(gca,'xtick',[],'ytick',[],'LineWidth',4,'fontsize',12,'FontWeight','bold');
set(gca,'Layer', 'top');

ylim([0 1])
hold on;

for ind = 1:n_old    
    figure(2)        
    plot(res_sigma_old(:,ind),'Color',color_ind,'LineWidth',4)
    hold on;
end
 
figure(2);
plot(mean_res_old,'Color',color_old,'LineWidth',6);

%% ........................................................................Plot figures for each individual
for i_young=1:size(Young, 2)
    
    CurrSubj = Young{i_young};
    
    sigma_individual = [];
    
    for i_cond = 1:size(Conditions,2)
        
        CurrCond = Conditions{i_cond};
        
        cd(Result_Dir);
        
        sigma_all = [CurrCond '_res_all_sigma_young'];
        load(sigma_all);
        
        sigma_condition = res_sigma_young(:,i_young);
        
        sigma_individual = [sigma_individual sigma_condition];
        
    end
    
    figure(1)
    
    set(gca,'xtick',[],'ytick',[],'LineWidth',4,'fontsize',12,'FontWeight','bold');
    
    ylim([0 1])
    hold on;
    
    plot(sigma_individual, 'LineWidth',4)
    
    cd(FigurePath);
    
    FigName = [CurrSubj '_rescale_sigma.png'];
    saveas(figure(1), FigName)
    
    close all
    
end

for i_old=1:size(Old, 2)
    
    CurrSubj = Old{i_old};
    
    sigma_individual = [];
    
    for i_cond = 1:size(Conditions,2)
        
        CurrCond = Conditions{i_cond};
        
        cd(Result_Dir);
        
        sigma_all = [CurrCond '_res_all_sigma_old'];
        load(sigma_all);
        
        sigma_condition = res_sigma_old(:,i_old);
        
        sigma_individual = [sigma_individual sigma_condition];
        
    end
    
    figure(1)
    
    set(gca,'xtick',[],'ytick',[],'LineWidth',4,'fontsize',12,'FontWeight','bold');
    
    ylim([0 1])
    hold on;
    
    plot(sigma_individual, 'LineWidth',4)
    
    cd(FigurePath);
    
    FigName = [CurrSubj '_rescale_sigma.png'];
    saveas(figure(1), FigName)
    
    close all
    
end