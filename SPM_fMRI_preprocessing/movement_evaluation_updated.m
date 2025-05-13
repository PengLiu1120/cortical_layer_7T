%% Script to plot and calculate the movementes in MoCo files.

%%% This script plots the movements displacements and the difference of movements between 
%%% the sample points of x,y,z axis from the MoCo files and saves the plot in a .png file. 

%%% written by Igor Tellez 31 July 2019
%%% contact: igor.tellez@ovgu.de

%%% updated by Juliane Doehler 14 August 2019
%%% contact: juliane.doehler@med.ovgu.de

%% clear matlab environment
clc; clear all; close all;

%% specifiy files and source folder
files = {'participant1','participant2'}

folder = '/home/doehlerj/Dokumente/PhD/sensemap/fMRI/data/'

%% iterate over all subjects
for fileIdx = 1:length(files)
    
    %% Load the output file (.mat) of the script "main01~" (DICOMtoNIFTI).
    foldername = sprintf('%s%s/resting_state/MotionParameters/', folder, files{fileIdx});
    meta = dir(fullfile(foldername,'*MoCoParam*.mat'));

    %% try to load file to evaluate and catch resulting load errors
    try
        temp = load(sprintf('%s%s',foldername,meta.name));
        success = 1;
    catch
        warning(sprintf('Cannot load dataset %s!',files{fileIdx}));
        success = 0;
    end

    %% if file is loaded procede
    if success == 1

        %% specifiy variables for evaluation plot
        fn = fieldnames(temp);
        movement_data = temp.(fn{1});
        id_val = extractBetween(folder, 'fMRI/', '/');
        idstr = string(id_val);
        plotname = strcat(idstr, '_Mov_plots');
        threshold = 2;    %value of threshold in milimeters


        %% Data processing 
        x = movement_data(:,1);
        y = movement_data(:,2);
        z = movement_data(:,3);
        meanX = mean(x);
        meanY = mean(y);
        meanZ = mean (z);
        x_axis = 1:length(x);
        x_ax = 1:length(x_axis)-1;

        xthreshold1 = meanX +threshold;
        xthreshold2 = meanX -threshold;
        aboveThresholdX1 = x>xthreshold1;
        aboveThresholdX2 = x<xthreshold2;

        ythreshold1 = meanY +threshold;
        ythreshold2 = meanY -threshold;
        aboveThresholdY1 = y>ythreshold1;
        aboveThresholdY2 = y<ythreshold2;

        zthreshold1 = meanZ +threshold;
        zthreshold2 = meanZ -threshold;
        aboveThresholdZ1 = z>zthreshold1;
        aboveThresholdZ2 = z<zthreshold2;

        % Difference between neighbors 

%         for i = 1:length(x_ax)
%             xdif(i) = x(i+1) - x(i); 
%             ydif(i) = y(i+1) - y(i); 
%             zdif(i) = z(i+1) - z(i);
%         end

        % Difference to the first image
        
         for i = 1:length(x_ax)
             xdif(i) = x(1) - x(i); 
             ydif(i) = y(1) - y(i); 
             zdif(i) = z(1) - z(i);
        end


        %% Plots

        meanXplot = ones(length(x),1)*meanX;
        meanYplot = ones(length(y),1)*meanY;
        meanZplot = ones(length(z),1)*meanZ;
        figure ('Name', 'Movements evaluation', 'NumberTitle', 'off')

        %X plot
        subplot (3,2,1) %define row, columns, plot number
        plot(x, 'LineWidth', 0.5)
        title('X-axis displacement')
        xlabel('Time')
        ylabel('Milimeters')
        hold all;
        plot(meanXplot, 'Color','k','LineStyle','--')
        %plot(x_axis(aboveThresholdX1),x(aboveThresholdX1),'r*')
        %plot(x_axis(aboveThresholdX2),x(aboveThresholdX2),'r*')
        xlim([1,length(x)])
        ylim([-4,4])
        %axis 'auto y'

        subplot(3,2,2)
        plot(x_ax, xdif)
        hold all;
        plot(x_ax(xdif>threshold), xdif(xdif>threshold), 'r*')
        title('\fontsize{8}X-axis difference between sample points')
        xlabel('Samples')
        ylabel('Milimeters')
        xlim([1,length(xdif)])
        ylim([-2,2])
        %axis 'auto y'
        set(gca,'FontSize',7)

        %Y plot
        subplot(3,2,3)
        plot(y)
        hold all;
        title('Y-axis displacement')
        xlabel('Time')
        ylabel('Milimeters')
        plot(meanYplot, 'Color', 'r', 'LineStyle', '--')
        %plot(x_axis(aboveThresholdY1),y(aboveThresholdY1),'r*')
        %plot(x_axis(aboveThresholdY2),y(aboveThresholdY2),'r*')
        xlim([1,length(x)])
        ylim([-4,4])
        %axis 'auto y'

        subplot(3,2,4)
        plot(x_ax, ydif)
        hold all;
        plot(x_ax(ydif>threshold), ydif(ydif>threshold), 'r*')
        title('\fontsize{8}Y-axis difference between sample points')
        xlabel('Samples')
        ylabel('Milimeters')
        xlim([1,length(ydif)])
        ylim([-2,2])
        %axis 'auto y'
        set(gca,'FontSize',7)

        %Z plot
        subplot(3,2,5)
        plot(z)
        title('Z-axis displacement')
        xlabel('Time')
        ylabel('Milimeters')
        hold all;
        plot(meanZplot, 'Color', 'r', 'LineStyle', '--')
        %plot(x_axis(aboveThresholdZ1),z(aboveThresholdZ1),'r*')
        %plot(x_axis(aboveThresholdZ2),z(aboveThresholdZ2),'r*')
        xlim([1,length(x)])
        ylim([-4,4])
        %axis 'auto y'

        subplot(3,2,6)
        plot(x_ax, zdif)
        hold all;
        plot(x_ax(zdif>threshold), zdif(zdif>threshold), 'r*')
        title('\fontsize{8}Z-axis difference between sample points')
        xlabel('Samples')
        ylabel('Milimeters')
        xlim([1,length(zdif)])
        ylim([-2,2])
        %axis 'auto y'
        set(gca,'FontSize',7)


        %% Save figure
        saveas(gca, sprintf('/home/doehlerj/Dokumente/PhD/sensemap/fMRI/resting_state/movement_evaluation/%s_MoCoParam_Evaluation_ylim.png',files{fileIdx}), 'png');
        disp('Done!')

    end

end
