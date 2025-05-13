%% Combine and percent signal change results
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
Result_Dir = '/media/pliu/LayerPRF/LayerMapping/08_LayerExtraction/Fourier';

FigurePath = '/home/pliu/Documents/LayerPRF/Figures/rescale_percent_signal_change_three_layers';

%% ........................................................................Set defaults
% .........................................................................Specify subjects
Young = {'frj712' 'gxo876' 'hby152' 'ijt563' 'kdy341' 'lpr469' 'nhm378' 'oms448' 'qet940' 'qxo538' 'unk742'};
Old = {'ajz367' 'bkn792' 'bmg520' 'cxc075' 'czg996' 'ggp057' 'gph998' 'iwq192' 'llh150' 'sst050'};

% .........................................................................Specify condtions, D2, D3 and D2+D3
Conditions = {'D2+D3' 'D2' 'D3'};

%% ........................................................................Specify group parameters
% .........................................................................Specify rescaled layer-specific percent signal change results
percent_signal_change_young='D3_percent_signal_change_young.mat';
load(percent_signal_change_young);
percent_signal_change_young = percent_signal_change_all;
percent_signal_change_old='D3_percent_signal_change_old.mat';
load(percent_signal_change_old);
percent_signal_change_old = percent_signal_change_all;

% .........................................................................Specify plot colours
colorspec_group = {'0000e1','e10000','868584'};

color_young = sscanf(colorspec_group{1},'%2x%2x%2x',[1 3])/255;
color_old = sscanf(colorspec_group{2},'%2x%2x%2x',[1 3])/255;
color_ind = sscanf(colorspec_group{3},'%2x%2x%2x',[1 3])/255;

% .........................................................................Specify group numbers
n_young = 11;
n_old = 10;

%% ........................................................................Calculate mean
mean_young = nanmean(percent_signal_change_young,2);
mean_old = nanmean(percent_signal_change_old,2);

%% ........................................................................Plot figures for young
figure(1)

set(gca,'xtick',[],'ytick',[],'LineWidth',4,'fontsize',12,'FontWeight','bold');
set(gca,'Layer', 'top');

xlim([0 4])
ylim([0 0.001])
hold on;

x_axis = [0.5, 2, 3.5];

for ind = 1:n_young    
    figure(1)        
    plot(x_axis,percent_signal_change_young(:,ind),'o','MarkerSize',6,'MarkerEdgeColor',color_ind,'MarkerFaceColor',color_ind)
    hold on;
end
 
figure(1);
plot(x_axis, mean_young,'o','MarkerSize',12,'MarkerEdgeColor',color_young,'MarkerFaceColor',color_young);

%% ........................................................................Plot figures for old
figure(2)

set(gca,'xtick',[],'ytick',[],'LineWidth',4,'fontsize',12,'FontWeight','bold');
set(gca,'Layer', 'top');

xlim([0 4])
ylim([0 0.001])
hold on;

x_axis = [0.5, 2, 3.5];

for ind = 1:n_old    
    figure(2)        
    plot(x_axis,percent_signal_change_old(:,ind),'o','MarkerSize',6,'MarkerEdgeColor',color_ind,'MarkerFaceColor',color_ind)
    hold on;
end
 
figure(2);
plot(x_axis,mean_old,'o','MarkerSize',12,'MarkerEdgeColor',color_old,'MarkerFaceColor',color_old);

%% ........................................................................Plot figures for each individual
for i_young=1:size(Young, 2)
    
    CurrSubj = Young{i_young};
    
    percent_signal_change_individual = [];
    
    for i_cond = 1:size(Conditions,2)
        
        CurrCond = Conditions{i_cond};
        
        cd(Result_Dir);
        
        percent_signal_change_all = [CurrCond '_res_percent_signal_change_young'];
        load(percent_signal_change_all);
        
        percent_signal_change_condition = res_percent_signal_change_young(:,i_young);
        
        percent_signal_change_individual = [percent_signal_change_individual percent_signal_change_condition];
        
    end
    
    figure(1)
    
    set(gca,'xtick',[],'ytick',[],'LineWidth',4,'fontsize',12,'FontWeight','bold');
    
    xlim([0 4])
    ylim([0 1])
    hold on;
    
    x_axis = [0.5, 2, 3.5];
    
    plot(x_axis,percent_signal_change_individual, '*','MarkerSize',6)
    
    cd(FigurePath);
    
    FigName = [CurrSubj '_rescale_percent_signal_change.png'];
    saveas(figure(1), FigName)
    
    close all
    
end

for i_old=1:size(Old, 2)
    
    CurrSubj = Old{i_old};
    
    percent_signal_change_individual = [];
    
    for i_cond = 1:size(Conditions,2)
        
        CurrCond = Conditions{i_cond};
        
        cd(Result_Dir);
        
        percent_signal_change_all = [CurrCond '_res_percent_signal_change_old'];
        load(percent_signal_change_all);
        
        percent_signal_change_condition = res_percent_signal_change_old(:,i_old);
        
        percent_signal_change_individual = [percent_signal_change_individual percent_signal_change_condition];
        
    end
    
    figure(1)
    
    set(gca,'xtick',[],'ytick',[],'LineWidth',4,'fontsize',12,'FontWeight','bold');
    
    xlim([0 4])
    ylim([0 1])
    hold on;
    
    x_axis = [0.5, 2, 3.5];
    
    plot(x_axis,percent_signal_change_individual, '*','MarkerSize',6)
    
    cd(FigurePath);
    
    FigName = [CurrSubj '_rescale_percent_signal_change.png'];
    saveas(figure(1), FigName)
    
    close all
    
end