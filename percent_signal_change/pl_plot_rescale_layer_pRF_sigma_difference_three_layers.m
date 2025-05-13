%% Plot rescaled layer-specific percent signal change difference results
% .........................................................................
% This script plots the rescaled layer-specific percent signal change
% difference results for young and old group respectively, and for each
% individual.
% .........................................................................
% Written by P.Liu
% Optimized by P.Liu
% Email: peng.liu@med.ovgu.de
% Last modified by P.Liu 14 May 2023
%% ........................................................................Tidy up
clear all
close all
clc

%% ........................................................................Set directory
Result_Dir = '/home/pliu/Documents/LayerPRF/Layer_Specific/08_LayerExtraction/pRF';
cd(Result_Dir);

FigurePath = '/home/pliu/Documents/LayerPRF/Figures/rescale_sigma_individual';

%% ........................................................................Set defaults
% .........................................................................Specify subjects
Young = {'frj712' 'gxo876' 'hby152' 'ijt563' 'kdy341' 'lpr469' 'nhm378' 'oms448' 'qet940' 'qxo538' 'unk742'};
Old = {'ajz367' 'bkn792' 'bmg520' 'cxc075' 'czg996' 'ggp057' 'gph998' 'iwq192' 'llh150' 'sst050'};

% .........................................................................Specify rescaled layer-specific percent signal change results
sigma_difference_young='sigma_difference_young.mat';
load(sigma_difference_young);
sigma_difference_old='sigma_difference_old.mat';
load(sigma_difference_old);

% .........................................................................Specify plot colours
colorspec_group = {'0000e1','e10000','868584'};

color_young = sscanf(colorspec_group{1},'%2x%2x%2x',[1 3])/255;
color_old = sscanf(colorspec_group{2},'%2x%2x%2x',[1 3])/255;
color_ind = sscanf(colorspec_group{3},'%2x%2x%2x',[1 3])/255;

% .........................................................................Specify group numbers
n_young = 11;
n_old = 10;

%% ........................................................................Calculate mean
mean_res_young = nanmean(sigma_difference_young,2);
mean_res_old = nanmean(sigma_difference_old,2);

%% ........................................................................Plot figures for young
figure(1)

set(gca,'xtick',[],'ytick',[],'LineWidth',4,'fontsize',12,'FontWeight','bold');
set(gca,'Layer', 'top');

xlim([0 4])
ylim([0 1])
hold on;

x_axis = [0.5, 2, 3.5];

for ind = 1:n_young
    figure(1)
    plot(x_axis,sigma_difference_young(:,ind),'o','MarkerSize',6,'MarkerEdgeColor',color_ind,'MarkerFaceColor',color_ind)
    hold on;
end

figure(1);
plot(x_axis, mean_res_young,'o','MarkerSize',12,'MarkerEdgeColor',color_young,'MarkerFaceColor',color_young);

%% ........................................................................Plot figures for old
figure(2)

set(gca,'xtick',[],'ytick',[],'LineWidth',4,'fontsize',12,'FontWeight','bold');
set(gca,'Layer', 'top');

xlim([0 4])
ylim([0 1])
hold on;

x_axis = [0.5, 2, 3.5];

for ind = 1:n_old
    figure(2)
    plot(x_axis,sigma_difference_old(:,ind),'o','MarkerSize',6,'MarkerEdgeColor',color_ind,'MarkerFaceColor',color_ind)
    hold on;
end

figure(2);
plot(x_axis, mean_res_old,'o','MarkerSize',12,'MarkerEdgeColor',color_old,'MarkerFaceColor',color_old);

%% ........................................................................Plot figures for each individual
for i_young=1:size(Young, 2)
    
    CurrSubj = Young{i_young};
    
    cd(Result_Dir);
    
    sigma_all = 'res_all_sigma_difference_young';
    load(sigma_all);
    
    sigma_individual = res_sigma_difference_young(:,i_young);
    
    figure(1)
    
    set(gca,'xtick',[],'ytick',[],'LineWidth',4,'fontsize',12,'FontWeight','bold');
    
    ylim([0 1])
    hold on;
    
    plot(sigma_individual,'b','LineWidth',4)
    
    cd(FigurePath);
    
    FigName = [CurrSubj '_rescale_sigma_difference.png'];
    saveas(figure(1), FigName)
    
    close all
    
end

for i_old=1:size(Old, 2)
    
    CurrSubj = Old{i_old};
    
    cd(Result_Dir);
    
    sigma_all = 'res_all_sigma_difference_old';
    load(sigma_all);
    
    sigma_individual = res_sigma_difference_old(:,i_old);
    
    figure(1)
    
    set(gca,'xtick',[],'ytick',[],'LineWidth',4,'fontsize',12,'FontWeight','bold');
    
    ylim([0 1])
    hold on;
    
    plot(sigma_individual,'r','LineWidth',4)
    
    cd(FigurePath);
    
    FigName = [CurrSubj '_rescale_sigma_difference.png'];
    saveas(figure(1), FigName)
    
    close all
    
end