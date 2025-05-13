%% Scatter plot for pRF mapping results
% This script plots scatter plot between Mu (vertical) and Sigma, Beta and
% Noise ceiling
% Written by P.Liu
% Last updated 12 Apr 2022 by P.Liu
%% ........................................................................Tidy up
clear all
close all
clc

%% ........................................................................Prepare
% .........................................................................Specify modelling path
RootPath = '/Users/pliu/Documents/LayerPRF/PRFmodelling';

% .........................................................................Specify subjects
Subjects = {'frj712' 'lpr469' 'qet940' 'qxo538' 'oms448'};
%'frj712' 'lpr469' 'qet940' 'qxo538' 'oms448'

% .........................................................................Specify condtions, D2, D3 and D2+D3
Conditions = {'D2' 'D3' 'D2+D3'};

% .........................................................................Specify ROIs, 3b hand and BA3b
ROI = {'3b_hand' '3b'};

% .........................................................................Specify fields
FldPRF = 'prf_';
FldFig = 'figures';

%% ........................................................................Scatter plot
% .........................................................................loop

for i_sub=1:size(Subjects, 2)
    
    CurrSubj = Subjects{i_sub};
    
    % .....................................................................Specify subject directory
    CurrSubjPath = fullfile(RootPath,CurrSubj);
    
    % .....................................................................Figure directory
    CurrFigPath = fullfile(RootPath,CurrSubj,FldFig);

    for i_cond = 1:size(Conditions,2)
        
        CurrCond = Conditions{i_cond};
        
        % .....................................................................Specify subject prf directory
        CurrPathPRF = fullfile(CurrSubjPath,[FldPRF CurrCond]);
        
        for i_ROI = 1:size(ROI,2)
            
            CurrROI = ROI{i_ROI};
            
            cd(CurrPathPRF)
            
            Map = ['lh_' CurrROI '_' CurrCond '_model.mat'];
            
            Vars = {'Values' 'TrueData'};
            
            Data = load(Map, Vars{:});
            
            nR2 = Data.TrueData(1,:);
            Mu = Data.TrueData(2,:);
            Sigma = Data.TrueData(3,:);
            Beta = Data.TrueData(4,:);
            NoiseCeiling = Data.TrueData(6,:);
            R2 = Data.TrueData(7,:);
                
            subplot(2,4,1)
            scatter(Mu,Sigma);
            xlabel('Mu','LineWidth',4,'FontSize',14);
            xlim([-1.5 1.5]);
            xticks(-1.5:1.5:1.5);
            ylabel('Sigma','LineWidth',4,'FontSize',14);
            ylim([0 3]);
            yticks(0:3:3);
            title('Mu vs Sigma','FontSize',14);
            
            subplot(2,4,2)
            scatter(Mu,Beta);
            xlabel('Mu','LineWidth',4,'FontSize',14);
            xlim([-1.5 1.5]);
            xticks(-1.5:1.5:1.5);
            ylabel('Beta','LineWidth',4,'FontSize',14);
            ylim([-0.1 0.1]);
            yticks(-0.1:0.1:0.1);
            title('Mu vs Beta','FontSize',14);
            
            subplot(2,4,3)
            scatter(Mu,NoiseCeiling);
            xlabel('Mu','LineWidth',4,'FontSize',14);
            xlim([-1.5 1.5]);
            xticks(-1.5:1.5:1.5);
            ylabel('NoiseCeiling','LineWidth',4,'FontSize',14);
            ylim([0 0.8]);
            yticks(0:0.8:0.8);
            title('Mu vs NoiseCeiling','FontSize',14);
            
            subplot(2,4,4)
            scatter(Sigma,R2);
            xlabel('Sigma','LineWidth',4,'FontSize',14);
            xlim([0 3]);
            xticks(0:3:3);
            ylabel('R^2','LineWidth',4,'FontSize',14);
            ylim([0 0.5]);
            yticks(0:0.5:0.5);
            title('Sigma vs R^2','FontSize',14);
            
            subplot(2,4,6)
            scatter(Sigma,Beta);
            xlabel('Sigma','LineWidth',4,'FontSize',14);
            xlim([0 3]);
            xticks(0:3:3);
            ylabel('Beta','LineWidth',4,'FontSize',14);
            ylim([-0.1 0.1]);
            yticks(-0.1:0.1:0.1);
            title('Sigma vs Beta','FontSize',14);
            
            subplot(2,4,7)
            scatter(Sigma,NoiseCeiling);
            xlabel('Sigma','LineWidth',4,'FontSize',14);
            xlim([0 3]);
            xticks(0:3:3);
            ylabel('NoiseCeiling','LineWidth',4,'FontSize',14);
            ylim([0 0.8]);
            yticks(0:0.8:0.8);
            title('Sigma vs NoiseCeiling','FontSize',14);
            
            subplot(2,4,8)
            scatter(Sigma,nR2);
            xlabel('Sigma','LineWidth',4,'FontSize',14);
            xlim([0 3]);
            xticks(0:3:3);
            ylabel('nR^2','LineWidth',4,'FontSize',14);
            ylim([0 1]);
            yticks(0:1:1);
            title('Sigma vs nR^2','FontSize',14);
            
            hold on
      
            set(gcf,'position',get(0, 'Screensize'));
            
            cd(CurrFigPath)
            
            FigName = ['Scatter_' CurrROI '_' CurrCond '.png'];
            
            saveas(gcf, FigName)
            
            close all
        
        end
        
    end
    
end