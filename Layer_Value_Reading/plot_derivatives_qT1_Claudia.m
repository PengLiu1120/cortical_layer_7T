
% this script reads in the h-file, binarizes the h-file, and multiplies the
% h-file with the T1-map, and reads out data

clear all;
close all;
clc;

%% General Specifications
subject = {'ajz367' 'bkn792' 'bmg520' 'cxc075' 'czg996' 'ggp057' 'gph998' 'iwq192' 'llh150' 'sst050'};
% Young: 'frj712' 'gxo876' 'hby152' 'ijt563' 'kdy341' 'lpr469' 'nhm378' 'oms448' 'qet940' 'qxo538' 'unk742'
% Old: 'ajz367' 'bkn792' 'bmg520' 'cxc075' 'czg996' 'ggp057' 'gph998' 'iwq192' 'llh150' 'sst050'


% definitions
% for laptop
% result_dir = '/Volumes/IKND/AG_Kuehn/Peng/LayerPRF/LayerMapping/06_Derivative';
T1_dir = 'C:\Users\Peng Liu admin\Documents\Nature_Neuroscience\Layer_Specific';

fileON="allsubjects_T1_old_hand_short.mat";
%layer_dir = '/media/doehlerj/WIP-B1/Documents/FINAL/T1_profiles/S1_Paper/results/curve_parameters'


colorspec_group = {'4D94FF','FF544D','dddddd','868584','FF7766'};
colorspec_layer = {'F4D7D7'};

color_grey = sscanf(colorspec_group{2},'%2x%2x%2x',[1 3])/255;
color_gy = sscanf(colorspec_group{1},'%2x%2x%2x',[1 3])/255;
color_gyl = sscanf(colorspec_group{4},'%2x%2x%2x',[1 3])/255;

%color_layerIV = sscanf(colorspec_layer{1},'%2x%2x%2x',[1 3])/255;

n_layers = 21; % specifies number of layers used in final plots
start_layer = 1;
end_layer = 21;
n_fingers = 5;

derivatives = [];

% Calculate mean deriv qT1 profile for whole group
mean_diff_S1_y=zeros(n_fingers,n_layers-1); % create empty variable to fill

load(sprintf('%s/%s',T1_dir,fileON))

table_allsubjects_ON = table_allsubjects_hand_short;

[num_rows, num_cols] = size(table_allsubjects_ON);

%load(sprintf('%s/idx_innerLayer.mat',layer_dir));
%load(sprintf('%s/idx_outerLayer.mat',layer_dir));
%load(sprintf('%s/idx_middleLayer.mat',layer_dir));

%inner=ind_inner
%outer=ind_outer
%middle=ind_middle

for n=1:num_rows
    diff_T1_ON = abs(gradient(table_allsubjects_ON(n,1:n_layers)));
    derivatives_ON(n,1:n_layers) = diff_T1_ON;
end

%% calculate mean derivatives
for deriv=1:n_layers % loop through all derivatives
    mean_diff_S1_y_ON(1,deriv)=nanmean(derivatives_ON(:,deriv));
end

figure(1)
%set(gca,'Xdir', 'reverse') % Sets x-axis to show superficial --> deep
% title('d2 + d3') % Sets figure title for the ROI
% xlabel({'WM <-- cortical depth --> pial'});
% ylabel('difference in qT1 [ms]');
set(gca, 'xtick', [1:4:21], 'ytick',[], 'LineWidth',4,'fontsize',28,'FontWeight','bold');
set(gca, 'Layer', 'top');
set(gca, 'View', [90 -90])
xlim([1 21])
ylim([0 120]) % Resets y-axis dimensions
hold on;

%patch('Faces',[1 2 3 4],'Vertices',[outer(7,2)+0.1 max(ylim); 21 max(ylim); 21 min(ylim); outer(7,2)+0.1 min(ylim)],'EdgeColor','none','FaceColor',color_grey);
%patch('Faces',[1 2 3 4],'Vertices',[middle(7,1)-0.1 max(ylim); middle(7,2)+0.1 max(ylim); middle(7,2)+0.1 min(ylim); middle(7,1)-0.1  min(ylim)],'EdgeColor','none','FaceColor',color_layerIV);
%patch('Faces',[1 2 3 4],'Vertices',[1 max(ylim); inner(7,1)-0.1 max(ylim); inner(7,1)-0.1 min(ylim); 1 min(ylim)],'EdgeColor','none','FaceColor',color_grey);


for ind = 1:num_rows    
    figure(1)        
    plot(derivatives_ON(ind,2:end),'Color',color_gyl,'LineWidth',6) % ensures that only current subject data is plotted 
    hold on;
end
 
figure(1);
% plot(mean_diff_S1_y_ON(1,:),'Color',color_gy,'LineWidth',12);
plot(mean_diff_S1_y_ON(1,:),'Color',color_grey,'LineWidth',12);