%% Plot rescaled layer-specific pRF size (sigma) results
% .........................................................................
% This script plots the rescaled layer-specific pRF size (sigma) results
% for young and old group respectively, and for each individual.
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
Result_Dir = 'C:\Users\Peng Liu admin\Documents\Nature_Neuroscience\Layer_Specific\08_LayerExtraction\pRF_Layer';
cd(Result_Dir);

FigureAllPath = 'C:\Users\Peng Liu admin\Documents\Nature_Neuroscience\Layer_Specific\08_LayerExtraction\pRF_Layer';
FigurePath = 'C:\Users\Peng Liu admin\Documents\Nature_Neuroscience\Layer_Specific\08_LayerExtraction\pRF_Layer';

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

%% ........................................................................
% .........................................................................Specify plot colours
colorspec_group = {'4D94FF','FF544D','868584'};

color_young = sscanf(colorspec_group{1},'%2x%2x%2x',[1 3])/255;
color_old = sscanf(colorspec_group{2},'%2x%2x%2x',[1 3])/255;
color_ind = sscanf(colorspec_group{3},'%2x%2x%2x',[1 3])/255;

% .........................................................................Specify group numbers
n_young = 11;
n_old = 10;

for i_cond = 1:size(Conditions,2)
    
    CurrCond = Conditions{i_cond};
    
    cd(Result_Dir);
    
    sigma_all = [CurrCond '_sigma_young'];
    load(sigma_all);
    
    sigma_all = [CurrCond '_sigma_old'];
    load(sigma_all);
    
    % .....................................................................Calculate mean
    mean_sigma_young = nanmean(sigma_young,2);
    mean_sigma_old = nanmean(sigma_old,2);
    
    %......................................................................Plot figures for young
    figure(1)
    
    set(gca,'xtick',[],'ytick',[],'LineWidth',4,'fontsize',12,'FontWeight','bold');
    set(gca,'Layer', 'top');
    
    ylim([0 1])
    hold on;
    
    for ind = 1:n_young
        figure(1)
        plot(sigma_young(:,ind),'Color',color_ind,'LineWidth',4)
        hold on;
    end
    
    figure(1);
    plot(mean_sigma_young,'Color',color_young,'LineWidth',6);
    
    cd(FigureAllPath);
    
    FigName = [CurrCond '_sigma_young.png'];
    saveas(figure(1), FigName)
    
    close all
    
    % .....................................................................Plot figures for old
    figure(2)
    
    set(gca,'xtick',[],'ytick',[],'LineWidth',4,'fontsize',12,'FontWeight','bold');
    set(gca,'Layer', 'top');
    
    ylim([0 1])
    hold on;
    
    for ind = 1:n_old
        figure(2)
        plot(sigma_old(:,ind),'Color',color_ind,'LineWidth',4)
        hold on;
    end
    
    figure(2);
    plot(mean_sigma_old,'Color',color_old,'LineWidth',6);
    
    cd(FigureAllPath);
    
    FigName = [CurrCond '_sigma_old.png'];
    saveas(figure(2), FigName)
    
    close all
    
end

%% ........................................................................Plot figures for each individual
for i_young=1:size(Young, 2)
    
    CurrSubj = Young{i_young};
    
    sigma_individual = [];
    
    for i_cond = 1:size(Conditions,2)
        
        CurrCond = Conditions{i_cond};
        
        cd(Result_Dir);
        
        sigma_all = [CurrCond '_sigma_young'];
        load(sigma_all);
        
        sigma_condition = sigma_young(:,i_young);
        
        sigma_individual = [sigma_individual sigma_condition];
        
    end
    
    figure(1)
    
    set(gca,'xtick',[],'ytick',[],'LineWidth',4,'fontsize',12,'FontWeight','bold');
    
    ylim([0 1])
    hold on;
    
    plot(sigma_individual, 'LineWidth',4)
    
    cd(FigurePath);
    
    FigName = [CurrSubj '_sigma.png'];
    saveas(figure(1), FigName)
    
    close all
    
end

for i_old=1:size(Old, 2)
    
    CurrSubj = Old{i_old};
    
    sigma_individual = [];
    
    for i_cond = 1:size(Conditions,2)
        
        CurrCond = Conditions{i_cond};
        
        cd(Result_Dir);
        
        sigma_all = [CurrCond '_sigma_old'];
        load(sigma_all);
        
        sigma_condition = sigma_old(:,i_old);
        
        sigma_individual = [sigma_individual sigma_condition];
        
    end
    
    figure(1)
    
    set(gca,'xtick',[],'ytick',[],'LineWidth',4,'fontsize',12,'FontWeight','bold');
    
    ylim([0 1])
    hold on;
    
    plot(sigma_individual, 'LineWidth',4)
    
    cd(FigurePath);
    
    FigName = [CurrSubj '_sigma.png'];
    saveas(figure(1), FigName)
    
    close all
    
end