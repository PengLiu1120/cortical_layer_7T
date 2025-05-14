% this script reads in the finger maps, binarizes the finger maps, and multiplies the
% finger maps with the T1-map, and reads out data

clear;

%% General Specifications
%%%%%%%%%%%%%%%%%%%%%%%%%
% subject definitions
%%%%%%%%%%%%%%%%%%%%%
fileID = fopen('/home/juli/projects/onehander/data/paths/all_IDs.txt','r');
formatSpec = '%3s';

A = textscan(fileID,formatSpec)

subject = A{1,1}'

%subject_folder = {'00','01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','20','21','22','23','24','25','26','27','28','29','30','31','32','33','34','35','36','37','38','39','40'};

subject_code = [3 3]; % 1=young, 2=old, 3=onehander

%gender = [1 2 2 1 1 2 1 1 2 1 1 1 1 2 1 1 1 2 2 1 2 2 1 2 2 1 1 2 1 2 2 2 2 1 1 2 2 2 2 1]; % 1=female, 2=male

script_dir = '/home/juli/projects/onehander/scripts';
result_dir = '/home/juli/projects/onehander/results/ReadOutLayerData_QSM_HandFaceFoot';
trs_dir = '/home/juli/projects/onehander/results/ReadOutLayerData_QSM_HandFaceFoot/thrs__better_for_nom';

addpath(script_dir)

% specifiy number of finger maps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_fingers = 3;

% some variables
%%%%%%%%%%%%%%%%
ind_y = 1;
ind_o = 1;

ly = [0 0 0];
lo = [0 0 0];

ly2 = [0 0 0];
lo2 = [0 0 0];

table_allsubjects_fingers_long = [];
table_allsubjects_fingers_short = [];

T = [];

subjects_young = {'dummy'};
subjects_old = {'dummy'};

NaN_subj = 1;

NaN_data_abs = {'dummy'};
NaN_data_code = [];

finger_sums_y = [0 0 0 0 0 0 0];
finger_sums_o = [0 0 0 0 0 0 0];

tmp = [];

% set flag pairwise to true to exlude missing values pairwise (i.e. per
% condition), default is false to exlude missing values casewise (i.e.
% exclude whole case from plotting if NaN values occur in at least one
% condition)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pairwise = false;
%pairwise = true;

% color specification
%%%%%%%%%%%%%%%%%%%%%
colorspec_group = {'000000','ff0000','dddddd','bbbbbb','767676'};
colorspec_layer = {'F4D7D7'};
colorspec = {'808080','556b2f','a0522d','191970','8b0000','808000','008000','3cb371','2f4f4f','008080','4682b4','9acd32','00008b','32cd32','daa520','7f007f','b03060','d2b48c','ff4500','00ced1','ff8c00','ffd700','6a5acd','0000cd','7cfc00','ba55d3','dc143c','00bfff','a020f0','ff69b4','f08080','d8bfd8','ff1493','ff00ff','1e90ff','f0e68c','b0e0e6','98fb98','7fffd4','bdb76b'};

color_grey = sscanf(colorspec_group{3},'%2x%2x%2x',[1 3])/255;
color_gy = sscanf(colorspec_group{1},'%2x%2x%2x',[1 3])/255;
color_go = sscanf(colorspec_group{2},'%2x%2x%2x',[1 3])/255;
color_gyl = sscanf(colorspec_group{4},'%2x%2x%2x',[1 3])/255;
color_gol = sscanf(colorspec_group{5},'%2x%2x%2x',[1 3])/255;

% specifies number of layers (used in final plots)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_layers = 21; 
start_layer = 1;
end_layer = 21;

% load files containing indices of extracted data driven layers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load('/media/doehlerj/WIP-B1/Documents/FINAL/T1_profiles/S1_Paper/results/curve_parameters/idx_innerLayer_y.mat');
% load('/media/doehlerj/WIP-B1/Documents/FINAL/T1_profiles/S1_Paper/results/curve_parameters/idx_middleLayer_y.mat');
% load('/media/doehlerj/WIP-B1/Documents/FINAL/T1_profiles/S1_Paper/results/curve_parameters/idx_outerLayer_y.mat');

% load('/media/doehlerj/WIP-B1/Documents/FINAL/T1_profiles/S1_Paper/results/curve_parameters/idx_innerLayer_o.mat');
% load('/media/doehlerj/WIP-B1/Documents/FINAL/T1_profiles/S1_Paper/results/curve_parameters/idx_middleLayer_o.mat');
% load('/media/doehlerj/WIP-B1/Documents/FINAL/T1_profiles/S1_Paper/results/curve_parameters/idx_outerLayer_o.mat');
% 
% inner=ind_inner
% outer=ind_outer
% middle=ind_middle

%% figure settings
%%%%%%%%%%%%%%%%%%
figure(1);
set(gcf, 'Position', [100 100 2100 500]);

subplot(1,7,1);
%set(gca,'Xdir', 'reverse') % Sets x-axis to show superficial --> deep
%set(gca,'Ydir', 'reverse')
title('hand action map'); % Sets figure title for the current subject - edit later
xlabel({'WM <-- cortical depth --> pial'});
ylabel('pQSM [ms]');
set(gca, 'xtick', 1:2:21);
% set(gca,'YMinorTick','on');
% set(gca,'XMinorTick','on');
set(gca, 'Layer', 'top');
set(gca,'View',[90 -90])
xlim([1 21]);
ylim([0 0.03]); % Resets y-axis dimensions
hold on;
%patch('Faces',[1 2 3 4],'Vertices',[19 max(ylim); 21 max(ylim); 21 min(ylim); 19 min(ylim)],'EdgeColor','none','FaceColor',color_grey);
%patch('Faces',[1 2 3 4],'Vertices',[middle(7,1) max(ylim); middle(7,2) max(ylim); middle(7,2) min(ylim); middle(7,1)  min(ylim)],'EdgeColor','none','FaceColor',color_layerIV);
%patch('Faces',[1 2 3 4],'Vertices',[1 max(ylim); 2 max(ylim); 2 min(ylim); 1 min(ylim)],'EdgeColor','none','FaceColor',color_grey);
% hold on;

subplot(1,7,2);
%set(gca,'Xdir', 'reverse') % Sets x-axis to show superficial --> deep
%set(gca,'Ydir', 'reverse')
title('face action map'); % Sets figure title for the current subject - edit later
xlabel({'WM <-- cortical depth --> pial'});
ylabel('pQSM [ppm]');
set(gca, 'xtick', 1:2:21);
% set(gca,'YMinorTick','on');
% set(gca,'XMinorTick','on');
set(gca, 'Layer', 'top');
set(gca,'View',[90 -90])
xlim([1 21]);
ylim([0 0.03]); % Resets y-axis dimensions
hold on;
%patch('Faces',[1 2 3 4],'Vertices',[19 max(ylim); 21 max(ylim); 21 min(ylim); 19 min(ylim)],'EdgeColor','none','FaceColor',color_grey);
%patch('Faces',[1 2 3 4],'Vertices',[middle(7,1) max(ylim); middle(7,2) max(ylim); middle(7,2) min(ylim); middle(7,1)  min(ylim)],'EdgeColor','none','FaceColor',color_layerIV);
%patch('Faces',[1 2 3 4],'Vertices',[1 max(ylim); 2 max(ylim); 2 min(ylim); 1 min(ylim)],'EdgeColor','none','FaceColor',color_grey);
% hold on;

subplot(1,7,3);
%set(gca,'Xdir', 'reverse') % Sets x-axis to show superficial --> deep
%set(gca,'Ydir', 'reverse')
title('foot action map'); % Sets figure title for the current subject - edit later
xlabel({'WM <-- cortical depth --> pial'});
ylabel('pQSM [ppm]');
set(gca, 'xtick', 1:2:21);
% set(gca,'YMinorTick','on');
% set(gca,'XMinorTick','on');
set(gca, 'Layer', 'top');
set(gca,'View',[90 -90])
xlim([1 21]);
ylim([0 0.03]); % Resets y-axis dimensions
hold on;
%patch('Faces',[1 2 3 4],'Vertices',[19 max(ylim); 21 max(ylim); 21 min(ylim); 19 min(ylim)],'EdgeColor','none','FaceColor',color_grey);
%patch('Faces',[1 2 3 4],'Vertices',[middle(7,1) max(ylim); middle(7,2) max(ylim); middle(7,2) min(ylim); middle(7,1)  min(ylim)],'EdgeColor','none','FaceColor',color_layerIV);
%patch('Faces',[1 2 3 4],'Vertices',[1 max(ylim); 2 max(ylim); 2 min(ylim); 1 min(ylim)],'EdgeColor','none','FaceColor',color_grey);
% hold on;

subplot(1,7,4);
%set(gca,'Xdir', 'reverse') % Sets x-axis to show superficial --> deep
%set(gca,'Ydir', 'reverse')
title('hand - 1st deriv'); % Sets figure title for the current subject - edit later
xlabel({'WM <-- cortical depth --> pial'});
ylabel('T1 [ms]');
set(gca, 'xtick', 1:2:21);
% set(gca,'YMinorTick','on');
% set(gca,'XMinorTick','on');
set(gca, 'Layer', 'top');
set(gca,'View',[90 -90])
xlim([1 21]);
ylim([0 150]); % Resets y-axis dimensions
hold on;
%patch('Faces',[1 2 3 4],'Vertices',[19 max(ylim); 21 max(ylim); 21 min(ylim); 19 min(ylim)],'EdgeColor','none','FaceColor',color_grey);
%patch('Faces',[1 2 3 4],'Vertices',[middle(7,1) max(ylim); middle(7,2) max(ylim); middle(7,2) min(ylim); middle(7,1)  min(ylim)],'EdgeColor','none','FaceColor',color_layerIV);
%patch('Faces',[1 2 3 4],'Vertices',[1 max(ylim); 2 max(ylim); 2 min(ylim); 1 min(ylim)],'EdgeColor','none','FaceColor',color_grey);
% hold on;

subplot(1,7,5);
%set(gca,'Xdir', 'reverse') % Sets x-axis to show superficial --> deep
%set(gca,'Ydir', 'reverse')
title('face - 1st deriv'); % Sets figure title for the current subject - edit later
xlabel({'WM <-- cortical depth --> pial'});
ylabel('T1 [ms]');
set(gca, 'xtick', 1:2:21);
% set(gca,'YMinorTick','on');
% set(gca,'XMinorTick','on');
set(gca, 'Layer', 'top');
set(gca,'View',[90 -90])
xlim([1 21]);
ylim([0 150]); % Resets y-axis dimensions
hold on;
%patch('Faces',[1 2 3 4],'Vertices',[19 max(ylim); 21 max(ylim); 21 min(ylim); 19 min(ylim)],'EdgeColor','none','FaceColor',color_grey);
%patch('Faces',[1 2 3 4],'Vertices',[middle(7,1) max(ylim); middle(7,2) max(ylim); middle(7,2) min(ylim); middle(7,1)  min(ylim)],'EdgeColor','none','FaceColor',color_layerIV);
%patch('Faces',[1 2 3 4],'Vertices',[1 max(ylim); 2 max(ylim); 2 min(ylim); 1 min(ylim)],'EdgeColor','none','FaceColor',color_grey);
% hold on;

subplot(1,7,6);
%set(gca,'Xdir', 'reverse') % Sets x-axis to show superficial --> deep
%set(gca,'Ydir', 'reverse')
title('foot - 1st deriv'); % Sets figure title for the current subject - edit later
xlabel({'WM <-- cortical depth --> pial'});
ylabel('T1 [ms]');
set(gca, 'xtick', 1:2:21);
% set(gca,'YMinorTick','on');
% set(gca,'XMinorTick','on');
set(gca, 'Layer', 'top');
set(gca,'View',[90 -90])
xlim([1 21]);
ylim([0 150]); % Resets y-axis dimensions
hold on;
%patch('Faces',[1 2 3 4],'Vertices',[19 max(ylim); 21 max(ylim); 21 min(ylim); 19 min(ylim)],'EdgeColor','none','FaceColor',color_grey);
%patch('Faces',[1 2 3 4],'Vertices',[middle(7,1) max(ylim); middle(7,2) max(ylim); middle(7,2) min(ylim); middle(7,1)  min(ylim)],'EdgeColor','none','FaceColor',color_layerIV);
%patch('Faces',[1 2 3 4],'Vertices',[1 max(ylim); 2 max(ylim); 2 min(ylim); 1 min(ylim)],'EdgeColor','none','FaceColor',color_grey);
% hold on;

subplot(1,7,7);
%set(gca,'Xdir', 'reverse') % Sets x-axis to show superficial --> deep
%set(gca,'Ydir', 'reverse')
title('OFF center'); % Sets figure title for the current subject - edit later
xlabel({'WM <-- cortical depth --> pial'});
ylabel('T1 [ms]');
set(gca, 'xtick', 1:2:21);
% set(gca,'YMinorTick','on');
% set(gca,'XMinorTick','on');
set(gca, 'Layer', 'top');
set(gca,'View',[90 -90])
xlim([1 21]);
ylim([1200 2800]); % Resets y-axis dimensions
hold on;
%patch('Faces',[1 2 3 4],'Vertices',[19 max(ylim); 21 max(ylim); 21 min(ylim); 19 min(ylim)],'EdgeColor','none','FaceColor',color_grey);
%patch('Faces',[1 2 3 4],'Vertices',[middle(7,1) max(ylim); middle(7,2) max(ylim); middle(7,2) min(ylim); middle(7,1)  min(ylim)],'EdgeColor','none','FaceColor',color_layerIV);
%patch('Faces',[1 2 3 4],'Vertices',[1 max(ylim); 2 max(ylim); 2 min(ylim); 1 min(ylim)],'EdgeColor','none','FaceColor',color_grey);
% hold on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2);
set(gcf, 'Position', [100 100 2100 500]);

subplot(1,7,1);
%set(gca,'Xdir', 'reverse') % Sets x-axis to show superficial --> deep
%set(gca,'Ydir', 'reverse')
title('hand action map'); % Sets figure title for the current subject - edit later
xlabel({'WM <-- cortical depth --> pial'});
ylabel('pQSM [ppm]');
set(gca, 'xtick', 1:2:21);
%set(gca,'YMinorTick','on');
%set(gca,'XMinorTick','on');
set(gca, 'Layer', 'top');
set(gca,'View',[90 -90])
xlim([1 21]);
ylim([0 0.03]); % Resets y-axis dimensions
hold on;
%patch('Faces',[1 2 3 4],'Vertices',[19 max(ylim); 21 max(ylim); 21 min(ylim); 19 min(ylim)],'EdgeColor','none','FaceColor',color_grey);
%patch('Faces',[1 2 3 4],'Vertices',[middle(7,1) max(ylim); middle(7,2) max(ylim); middle(7,2) min(ylim); middle(7,1)  min(ylim)],'EdgeColor','none','FaceColor',color_layerIV);
%patch('Faces',[1 2 3 4],'Vertices',[1 max(ylim); 2 max(ylim); 2 min(ylim); 1 min(ylim)],'EdgeColor','none','FaceColor',color_grey);
% hold on;

subplot(1,7,2);
%set(gca,'Xdir', 'reverse') % Sets x-axis to show superficial --> deep
%set(gca,'Ydir', 'reverse')
title('face action map'); % Sets figure title for the current subject - edit later
xlabel({'WM <-- cortical depth --> pial'});
ylabel('pQSM [ppm]');
set(gca, 'xtick', 1:2:21);
% set(gca,'YMinorTick','on');
% set(gca,'XMinorTick','on');
set(gca,'View',[90 -90])
set(gca, 'Layer', 'top');
xlim([1 21]);
ylim([0 0.03]); % Resets y-axis dimensions
hold on;
%patch('Faces',[1 2 3 4],'Vertices',[19 max(ylim); 21 max(ylim); 21 min(ylim); 19 min(ylim)],'EdgeColor','none','FaceColor',color_grey);
%patch('Faces',[1 2 3 4],'Vertices',[middle(7,1) max(ylim); middle(7,2) max(ylim); middle(7,2) min(ylim); middle(7,1)  min(ylim)],'EdgeColor','none','FaceColor',color_layerIV);
%patch('Faces',[1 2 3 4],'Vertices',[1 max(ylim); 2 max(ylim); 2 min(ylim); 1 min(ylim)],'EdgeColor','none','FaceColor',color_grey);
% hold on;

subplot(1,7,3);
%set(gca,'Xdir', 'reverse') % Sets x-axis to show superficial --> deep
%set(gca,'Ydir', 'reverse')
title('foot action map'); % Sets figure title for the current subject - edit later
xlabel({'WM <-- cortical depth --> pial'});
ylabel('pQSM [ppm]');
set(gca, 'xtick', 1:2:21);
% set(gca,'YMinorTick','on');
% set(gca,'XMinorTick','on');
set(gca, 'Layer', 'top');
set(gca,'View',[90 -90])
xlim([1 21]);
ylim([0 0.03]); % Resets y-axis dimensions
hold on;
%patch('Faces',[1 2 3 4],'Vertices',[19 max(ylim); 21 max(ylim); 21 min(ylim); 19 min(ylim)],'EdgeColor','none','FaceColor',color_grey);
%patch('Faces',[1 2 3 4],'Vertices',[middle(7,1) max(ylim); middle(7,2) max(ylim); middle(7,2) min(ylim); middle(7,1)  min(ylim)],'EdgeColor','none','FaceColor',color_layerIV);
%patch('Faces',[1 2 3 4],'Vertices',[1 max(ylim); 2 max(ylim); 2 min(ylim); 1 min(ylim)],'EdgeColor','none','FaceColor',color_grey);
% hold on;

subplot(1,7,4);
%set(gca,'Xdir', 'reverse') % Sets x-axis to show superficial --> deep
%set(gca,'Ydir', 'reverse')
title('hand - 1st deriv'); % Sets figure title for the current subject - edit later
xlabel({'WM <-- cortical depth --> pial'});
ylabel('T1 [ms]');
set(gca, 'xtick', 1:2:21);
% set(gca,'YMinorTick','on');
% set(gca,'XMinorTick','on');
set(gca, 'Layer', 'top');
set(gca,'View',[90 -90])
xlim([1 21]);
ylim([0 150]); % Resets y-axis dimensions
hold on;
%patch('Faces',[1 2 3 4],'Vertices',[19 max(ylim); 21 max(ylim); 21 min(ylim); 19 min(ylim)],'EdgeColor','none','FaceColor',color_grey);
%patch('Faces',[1 2 3 4],'Vertices',[middle(7,1) max(ylim); middle(7,2) max(ylim); middle(7,2) min(ylim); middle(7,1)  min(ylim)],'EdgeColor','none','FaceColor',color_layerIV);
%patch('Faces',[1 2 3 4],'Vertices',[1 max(ylim); 2 max(ylim); 2 min(ylim); 1 min(ylim)],'EdgeColor','none','FaceColor',color_grey);
% hold on;

subplot(1,7,5);
%set(gca,'Xdir', 'reverse') % Sets x-axis to show superficial --> deep
%set(gca,'Ydir', 'reverse')
title('face - 1st deriv'); % Sets figure title for the current subject - edit later
xlabel({'WM <-- cortical depth --> pial'});
ylabel('T1 [ms]');
set(gca, 'xtick', 1:2:21);
% set(gca,'YMinorTick','on');
% set(gca,'XMinorTick','on');
set(gca, 'Layer', 'top');
set(gca,'View',[90 -90])
xlim([1 21]);
ylim([0 150]); % Resets y-axis dimensions
hold on;
%patch('Faces',[1 2 3 4],'Vertices',[19 max(ylim); 21 max(ylim); 21 min(ylim); 19 min(ylim)],'EdgeColor','none','FaceColor',color_grey);
%patch('Faces',[1 2 3 4],'Vertices',[middle(7,1) max(ylim); middle(7,2) max(ylim); middle(7,2) min(ylim); middle(7,1)  min(ylim)],'EdgeColor','none','FaceColor',color_layerIV);
%patch('Faces',[1 2 3 4],'Vertices',[1 max(ylim); 2 max(ylim); 2 min(ylim); 1 min(ylim)],'EdgeColor','none','FaceColor',color_grey);
% hold on;

subplot(1,7,6);
%set(gca,'Xdir', 'reverse') % Sets x-axis to show superficial --> deep
%set(gca,'Ydir', 'reverse')
title('foot - 1st deriv'); % Sets figure title for the current subject - edit later
xlabel({'WM <-- cortical depth --> pial'});
ylabel('T1 [ms]');
set(gca, 'xtick', 1:2:21);
% set(gca,'YMinorTick','on');
% set(gca,'XMinorTick','on');
set(gca, 'Layer', 'top');
set(gca,'View',[90 -90])
xlim([1 21]);
ylim([0 150]); % Resets y-axis dimensions
hold on;
%patch('Faces',[1 2 3 4],'Vertices',[19 max(ylim); 21 max(ylim); 21 min(ylim); 19 min(ylim)],'EdgeColor','none','FaceColor',color_grey);
%patch('Faces',[1 2 3 4],'Vertices',[middle(7,1) max(ylim); middle(7,2) max(ylim); middle(7,2) min(ylim); middle(7,1)  min(ylim)],'EdgeColor','none','FaceColor',color_layerIV);
%patch('Faces',[1 2 3 4],'Vertices',[1 max(ylim); 2 max(ylim); 2 min(ylim); 1 min(ylim)],'EdgeColor','none','FaceColor',color_grey);
% hold on;

subplot(1,7,7);
%set(gca,'Xdir', 'reverse') % Sets x-axis to show superficial --> deep
%set(gca,'Ydir', 'reverse')
title('OFF center'); % Sets figure title for the current subject - edit later
xlabel({'WM <-- cortical depth --> pial'});
ylabel('T1 [ms]');
set(gca, 'xtick', 1:2:21);
% set(gca,'YMinorTick','on');
% set(gca,'XMinorTick','on');
set(gca, 'Layer', 'top');
set(gca,'View',[90 -90])
xlim([1 21]);
ylim([1200 2800]); % Resets y-axis dimensions
hold on;
%patch('Faces',[1 2 3 4],'Vertices',[19 max(ylim); 21 max(ylim); 21 min(ylim); 19 min(ylim)],'EdgeColor','none','FaceColor',color_grey);
%patch('Faces',[1 2 3 4],'Vertices',[middle(7,1) max(ylim); middle(7,2) max(ylim); middle(7,2) min(ylim); middle(7,1)  min(ylim)],'EdgeColor','none','FaceColor',color_layerIV);
%patch('Faces',[1 2 3 4],'Vertices',[1 max(ylim); 2 max(ylim); 2 min(ylim); 1 min(ylim)],'EdgeColor','none','FaceColor',color_grey);
% hold on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(3);
set(gcf, 'Position', [100 100 2100 500]);

subplot(1,7,1);
%set(gca,'Xdir', 'reverse') % Sets x-axis to show superficial --> deep
set(gca,'Ydir', 'reverse')
title('hand action map'); % Sets figure title for the current subject - edit later
xlabel({'WM <-- cortical depth --> pial'});
ylabel('T1 [ms]');
set(gca, 'xtick', 1:2:21);
% set(gca,'YMinorTick','on');
% set(gca,'XMinorTick','on');
set(gca, 'Layer', 'top');
set(gca,'View',[90 -90])
xlim([1 21]);
ylim([1200 2800]); % Resets y-axis dimensions
hold on;
%patch('Faces',[1 2 3 4],'Vertices',[19 max(ylim); 21 max(ylim); 21 min(ylim); 19 min(ylim)],'EdgeColor','none','FaceColor',color_grey);
%patch('Faces',[1 2 3 4],'Vertices',[middle(7,1) max(ylim); middle(7,2) max(ylim); middle(7,2) min(ylim); middle(7,1)  min(ylim)],'EdgeColor','none','FaceColor',color_layerIV);
%patch('Faces',[1 2 3 4],'Vertices',[1 max(ylim); 2 max(ylim); 2 min(ylim); 1 min(ylim)],'EdgeColor','none','FaceColor',color_grey);
% hold on;

subplot(1,7,2);
%set(gca,'Xdir', 'reverse') % Sets x-axis to show superficial --> deep
set(gca,'Ydir', 'reverse')
title('face action map'); % Sets figure title for the current subject - edit later
xlabel({'WM <-- cortical depth --> pial'});
ylabel('T1 [ms]');
set(gca, 'xtick', 1:2:21);
% set(gca,'YMinorTick','on');
% set(gca,'XMinorTick','on');
set(gca, 'Layer', 'top');
set(gca,'View',[90 -90])
xlim([1 21]);
ylim([1200 2800]); % Resets y-axis dimensions
hold on;
%patch('Faces',[1 2 3 4],'Vertices',[19 max(ylim); 21 max(ylim); 21 min(ylim); 19 min(ylim)],'EdgeColor','none','FaceColor',color_grey);
%patch('Faces',[1 2 3 4],'Vertices',[middle(7,1) max(ylim); middle(7,2) max(ylim); middle(7,2) min(ylim); middle(7,1)  min(ylim)],'EdgeColor','none','FaceColor',color_layerIV);
%patch('Faces',[1 2 3 4],'Vertices',[1 max(ylim); 2 max(ylim); 2 min(ylim); 1 min(ylim)],'EdgeColor','none','FaceColor',color_grey);
% hold on;

subplot(1,7,3);
%set(gca,'Xdir', 'reverse') % Sets x-axis to show superficial --> deep
set(gca,'Ydir', 'reverse')
title('foot action map'); % Sets figure title for the current subject - edit later
xlabel({'WM <-- cortical depth --> pial'});
ylabel('T1 [ms]');
set(gca, 'xtick', 1:2:21);
% set(gca,'YMinorTick','on');
% set(gca,'XMinorTick','on');
set(gca, 'Layer', 'top');
set(gca,'View',[90 -90])
xlim([1 21]);
ylim([1200 2800]); % Resets y-axis dimensions
hold on;
%patch('Faces',[1 2 3 4],'Vertices',[19 max(ylim); 21 max(ylim); 21 min(ylim); 19 min(ylim)],'EdgeColor','none','FaceColor',color_grey);
%patch('Faces',[1 2 3 4],'Vertices',[middle(7,1) max(ylim); middle(7,2) max(ylim); middle(7,2) min(ylim); middle(7,1)  min(ylim)],'EdgeColor','none','FaceColor',color_layerIV);
%patch('Faces',[1 2 3 4],'Vertices',[1 max(ylim); 2 max(ylim); 2 min(ylim); 1 min(ylim)],'EdgeColor','none','FaceColor',color_grey);
% hold on;

subplot(1,7,4);
%set(gca,'Xdir', 'reverse') % Sets x-axis to show superficial --> deep
%set(gca,'Ydir', 'reverse')
title('hand - 1st deriv'); % Sets figure title for the current subject - edit later
xlabel({'WM <-- cortical depth --> pial'});
ylabel('T1 [ms]');
set(gca, 'xtick', 1:2:21);
% set(gca,'YMinorTick','on');
% set(gca,'XMinorTick','on');
set(gca, 'Layer', 'top');
set(gca,'View',[90 -90])
xlim([1 21]);
ylim([0 150]); % Resets y-axis dimensions
hold on;
%patch('Faces',[1 2 3 4],'Vertices',[19 max(ylim); 21 max(ylim); 21 min(ylim); 19 min(ylim)],'EdgeColor','none','FaceColor',color_grey);
%patch('Faces',[1 2 3 4],'Vertices',[middle(7,1) max(ylim); middle(7,2) max(ylim); middle(7,2) min(ylim); middle(7,1)  min(ylim)],'EdgeColor','none','FaceColor',color_layerIV);
%patch('Faces',[1 2 3 4],'Vertices',[1 max(ylim); 2 max(ylim); 2 min(ylim); 1 min(ylim)],'EdgeColor','none','FaceColor',color_grey);
% hold on;

subplot(1,7,5);
%set(gca,'Xdir', 'reverse') % Sets x-axis to show superficial --> deep
%set(gca,'Ydir', 'reverse')
title('face - 1st deriv'); % Sets figure title for the current subject - edit later
xlabel({'WM <-- cortical depth --> pial'});
ylabel('T1 [ms]');
set(gca, 'xtick', 1:2:21);
% set(gca,'YMinorTick','on');
% set(gca,'XMinorTick','on');
set(gca, 'Layer', 'top');
set(gca,'View',[90 -90])
xlim([1 21]);
ylim([0 150]); % Resets y-axis dimensions
hold on;
%patch('Faces',[1 2 3 4],'Vertices',[19 max(ylim); 21 max(ylim); 21 min(ylim); 19 min(ylim)],'EdgeColor','none','FaceColor',color_grey);
%patch('Faces',[1 2 3 4],'Vertices',[middle(7,1) max(ylim); middle(7,2) max(ylim); middle(7,2) min(ylim); middle(7,1)  min(ylim)],'EdgeColor','none','FaceColor',color_layerIV);
%patch('Faces',[1 2 3 4],'Vertices',[1 max(ylim); 2 max(ylim); 2 min(ylim); 1 min(ylim)],'EdgeColor','none','FaceColor',color_grey);
% hold on;

subplot(1,7,6);
%set(gca,'Xdir', 'reverse') % Sets x-axis to show superficial --> deep
%set(gca,'Ydir', 'reverse')
title('foot - 1st deriv'); % Sets figure title for the current subject - edit later
xlabel({'WM <-- cortical depth --> pial'});
ylabel('T1 [ms]');
set(gca, 'xtick', 1:2:21);
% set(gca,'YMinorTick','on');
% set(gca,'XMinorTick','on');
set(gca, 'Layer', 'top');
set(gca,'View',[90 -90])
xlim([1 21]);
ylim([0 150]); % Resets y-axis dimensions
hold on;
%patch('Faces',[1 2 3 4],'Vertices',[19 max(ylim); 21 max(ylim); 21 min(ylim); 19 min(ylim)],'EdgeColor','none','FaceColor',color_grey);
%patch('Faces',[1 2 3 4],'Vertices',[middle(7,1) max(ylim); middle(7,2) max(ylim); middle(7,2) min(ylim); middle(7,1)  min(ylim)],'EdgeColor','none','FaceColor',color_layerIV);
%patch('Faces',[1 2 3 4],'Vertices',[1 max(ylim); 2 max(ylim); 2 min(ylim); 1 min(ylim)],'EdgeColor','none','FaceColor',color_grey);
% hold on;

subplot(1,7,7);
%set(gca,'Xdir', 'reverse') % Sets x-axis to show superficial --> deep
set(gca,'Ydir', 'reverse')
title('OFF center'); % Sets figure title for the current subject - edit later
xlabel({'WM <-- cortical depth --> pial'});
ylabel('T1 [ms]');
set(gca, 'xtick', 1:2:21);
% set(gca,'YMinorTick','on');
% set(gca,'XMinorTick','on');
set(gca, 'Layer', 'top');
set(gca,'View',[90 -90])
xlim([1 21]);
ylim([1200 2800]); % Resets y-axis dimensions
hold on;
%patch('Faces',[1 2 3 4],'Vertices',[19 max(ylim); 21 max(ylim); 21 min(ylim); 19 min(ylim)],'EdgeColor','none','FaceColor',color_grey);
%patch('Faces',[1 2 3 4],'Vertices',[middle(7,1) max(ylim); middle(7,2) max(ylim); middle(7,2) min(ylim); middle(7,1)  min(ylim)],'EdgeColor','none','FaceColor',color_layerIV);
%patch('Faces',[1 2 3 4],'Vertices',[1 max(ylim); 2 max(ylim); 2 min(ylim); 1 min(ylim)],'EdgeColor','none','FaceColor',color_grey);
% hold on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(4);
set(gcf, 'Position', [100 100 2100 500]);

subplot(1,7,1);
%set(gca,'Xdir', 'reverse') % Sets x-axis to show superficial --> deep
set(gca,'Ydir', 'reverse')
title('hand action map'); % Sets figure title for the current subject - edit later
xlabel({'WM <-- cortical depth --> pial'});
ylabel('T1 [ms]');
set(gca, 'xtick', 1:2:21);
% set(gca,'YMinorTick','on');
% set(gca,'XMinorTick','on');
set(gca, 'Layer', 'top');
set(gca,'View',[90 -90])
xlim([1 21]);
ylim([1200 2800]); % Resets y-axis dimensions
hold on;
%patch('Faces',[1 2 3 4],'Vertices',[19 max(ylim); 21 max(ylim); 21 min(ylim); 19 min(ylim)],'EdgeColor','none','FaceColor',color_grey);
%patch('Faces',[1 2 3 4],'Vertices',[middle(7,1) max(ylim); middle(7,2) max(ylim); middle(7,2) min(ylim); middle(7,1)  min(ylim)],'EdgeColor','none','FaceColor',color_layerIV);
%patch('Faces',[1 2 3 4],'Vertices',[1 max(ylim); 2 max(ylim); 2 min(ylim); 1 min(ylim)],'EdgeColor','none','FaceColor',color_grey);
% hold on;

subplot(1,7,2);
%set(gca,'Xdir', 'reverse') % Sets x-axis to show superficial --> deep
set(gca,'Ydir', 'reverse')
title('face action map'); % Sets figure title for the current subject - edit later
xlabel({'WM <-- cortical depth --> pial'});
ylabel('T1 [ms]');
set(gca, 'xtick', 1:2:21);
% set(gca,'YMinorTick','on');
% set(gca,'XMinorTick','on');
set(gca, 'Layer', 'top');
set(gca,'View',[90 -90])
xlim([1 21]);
ylim([1200 2800]); % Resets y-axis dimensions
hold on;
%patch('Faces',[1 2 3 4],'Vertices',[19 max(ylim); 21 max(ylim); 21 min(ylim); 19 min(ylim)],'EdgeColor','none','FaceColor',color_grey);
%patch('Faces',[1 2 3 4],'Vertices',[middle(7,1) max(ylim); middle(7,2) max(ylim); middle(7,2) min(ylim); middle(7,1)  min(ylim)],'EdgeColor','none','FaceColor',color_layerIV);
%patch('Faces',[1 2 3 4],'Vertices',[1 max(ylim); 2 max(ylim); 2 min(ylim); 1 min(ylim)],'EdgeColor','none','FaceColor',color_grey);
% hold on;

subplot(1,7,3);
%set(gca,'Xdir', 'reverse') % Sets x-axis to show superficial --> deep
set(gca,'Ydir', 'reverse')
title('foot action map'); % Sets figure title for the current subject - edit later
xlabel({'WM <-- cortical depth --> pial'});
ylabel('T1 [ms]');
set(gca, 'xtick', 1:2:21);
% set(gca,'YMinorTick','on');
% set(gca,'XMinorTick','on');
set(gca, 'Layer', 'top');
set(gca,'View',[90 -90])
xlim([1 21]);
ylim([1200 2800]); % Resets y-axis dimensions
hold on;
%patch('Faces',[1 2 3 4],'Vertices',[19 max(ylim); 21 max(ylim); 21 min(ylim); 19 min(ylim)],'EdgeColor','none','FaceColor',color_grey);
%patch('Faces',[1 2 3 4],'Vertices',[middle(7,1) max(ylim); middle(7,2) max(ylim); middle(7,2) min(ylim); middle(7,1)  min(ylim)],'EdgeColor','none','FaceColor',color_layerIV);
%patch('Faces',[1 2 3 4],'Vertices',[1 max(ylim); 2 max(ylim); 2 min(ylim); 1 min(ylim)],'EdgeColor','none','FaceColor',color_grey);
% hold on;

subplot(1,7,4);
%set(gca,'Xdir', 'reverse') % Sets x-axis to show superficial --> deep
%set(gca,'Ydir', 'reverse')
title('hand - 1st deriv'); % Sets figure title for the current subject - edit later
xlabel({'WM <-- cortical depth --> pial'});
ylabel('T1 [ms]');
set(gca, 'xtick', 1:2:21);
% set(gca,'YMinorTick','on');
% set(gca,'XMinorTick','on');
set(gca, 'Layer', 'top');
set(gca,'View',[90 -90])
xlim([1 21]);
ylim([0 150]); % Resets y-axis dimensions
hold on;
%patch('Faces',[1 2 3 4],'Vertices',[19 max(ylim); 21 max(ylim); 21 min(ylim); 19 min(ylim)],'EdgeColor','none','FaceColor',color_grey);
%patch('Faces',[1 2 3 4],'Vertices',[middle(7,1) max(ylim); middle(7,2) max(ylim); middle(7,2) min(ylim); middle(7,1)  min(ylim)],'EdgeColor','none','FaceColor',color_layerIV);
%patch('Faces',[1 2 3 4],'Vertices',[1 max(ylim); 2 max(ylim); 2 min(ylim); 1 min(ylim)],'EdgeColor','none','FaceColor',color_grey);
% hold on;

subplot(1,7,5);
%set(gca,'Xdir', 'reverse') % Sets x-axis to show superficial --> deep
%set(gca,'Ydir', 'reverse')
title('face - 1st deriv'); % Sets figure title for the current subject - edit later
xlabel({'WM <-- cortical depth --> pial'});
ylabel('T1 [ms]');
set(gca, 'xtick', 1:2:21);
% set(gca,'YMinorTick','on');
% set(gca,'XMinorTick','on');
set(gca, 'Layer', 'top');
set(gca,'View',[90 -90])
xlim([1 21]);
ylim([0 150]); % Resets y-axis dimensions
hold on;
%patch('Faces',[1 2 3 4],'Vertices',[19 max(ylim); 21 max(ylim); 21 min(ylim); 19 min(ylim)],'EdgeColor','none','FaceColor',color_grey);
%patch('Faces',[1 2 3 4],'Vertices',[middle(7,1) max(ylim); middle(7,2) max(ylim); middle(7,2) min(ylim); middle(7,1)  min(ylim)],'EdgeColor','none','FaceColor',color_layerIV);
%patch('Faces',[1 2 3 4],'Vertices',[1 max(ylim); 2 max(ylim); 2 min(ylim); 1 min(ylim)],'EdgeColor','none','FaceColor',color_grey);
% hold on;

subplot(1,7,6);
%set(gca,'Xdir', 'reverse') % Sets x-axis to show superficial --> deep
%set(gca,'Ydir', 'reverse')
title('foot - 1st deriv'); % Sets figure title for the current subject - edit later
xlabel({'WM <-- cortical depth --> pial'});
ylabel('T1 [ms]');
set(gca, 'xtick', 1:2:21);
% set(gca,'YMinorTick','on');
% set(gca,'XMinorTick','on');
set(gca, 'Layer', 'top');
set(gca,'View',[90 -90])
xlim([1 21]);
ylim([0 150]); % Resets y-axis dimensions
hold on;
%patch('Faces',[1 2 3 4],'Vertices',[19 max(ylim); 21 max(ylim); 21 min(ylim); 19 min(ylim)],'EdgeColor','none','FaceColor',color_grey);
%patch('Faces',[1 2 3 4],'Vertices',[middle(7,1) max(ylim); middle(7,2) max(ylim); middle(7,2) min(ylim); middle(7,1)  min(ylim)],'EdgeColor','none','FaceColor',color_layerIV);
%patch('Faces',[1 2 3 4],'Vertices',[1 max(ylim); 2 max(ylim); 2 min(ylim); 1 min(ylim)],'EdgeColor','none','FaceColor',color_grey);
% hold on;

subplot(1,7,7);
%set(gca,'Xdir', 'reverse') % Sets x-axis to show superficial --> deep
set(gca,'Ydir', 'reverse')
title('OFF center'); % Sets figure title for the current subject - edit later
xlabel({'WM <-- cortical depth --> pial'});
ylabel('T1 [ms]');
set(gca, 'xtick', 1:2:21);
% set(gca,'YMinorTick','on');
% set(gca,'XMinorTick','on');
set(gca, 'Layer', 'top');
set(gca,'View',[90 -90])
xlim([1 21]);
ylim([1200 2800]); % Resets y-axis dimensions
hold on;
%patch('Faces',[1 2 3 4],'Vertices',[19 max(ylim); 21 max(ylim); 21 min(ylim); 19 min(ylim)],'EdgeColor','none','FaceColor',color_grey);
%patch('Faces',[1 2 3 4],'Vertices',[middle(7,1) max(ylim); middle(7,2) max(ylim); middle(7,2) min(ylim); middle(7,1)  min(ylim)],'EdgeColor','none','FaceColor',color_layerIV);
%patch('Faces',[1 2 3 4],'Vertices',[1 max(ylim); 2 max(ylim); 2 min(ylim); 1 min(ylim)],'EdgeColor','none','FaceColor',color_grey);
% hold on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hemi = {'left','right'};

%res = {'orig','resampled'};
res = {'orig'};

%% start loop to: extract T1 values and finger maps for all subjects (one iteration per subject subject), store and plot data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ind = 1:length(subject)

    for h = 1:length(hemi)

        for r = 1:length(res)

             tmp = [];

             if ind == 1
                cd(result_dir)
             else
                cd(result_dir);
             end

             try
                load(sprintf('%s_pQSM_mface_21layer_mean_thrs_%s_%s_70_minROI.mat',subject{ind},hemi{h},res{r}));
                load(sprintf('%s_pQSM_mhand_21layer_mean_thrs_%s_%s_70_minROI.mat',subject{ind},hemi{h},res{r}));
                load(sprintf('%s_pQSM_mfoot_21layer_mean_thrs_%s_%s_70_minROI.mat',subject{ind},hemi{h},res{r}));
                actionmaps = 1;
             catch
                warning('No action maps available. Trying next subject');
                actionmaps = 0;
             end
        
             if actionmaps == 1
                 % add all data to sample table
                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                 table_allsubjects_fingers_long = [table_allsubjects_fingers_long; mean_T1_D1_21;mean_T1_D2_21;mean_T1_D3_21];
                 %table_allsubjects_fingers_short(ind,:) = [transpose(mean_T1_D1(:,1)) transpose(mean_T1_D2(:,1)) transpose(mean_T1_D3(:,1)) transpose(mean_T1_D4(:,1)) transpose(mean_T1_D5(:,1)) transpose(OFFcenter_mean) subject_code(ind) ind transpose(mean_T1_D1(1,6)) transpose(mean_T1_D2(1,6)) transpose(mean_T1_D3(1,6)) transpose(mean_T1_D4(1,6)) transpose(mean_T1_D5(1,6)) gender(ind)];
                 table_allsubjects_fingers_short(ind,:) = [transpose(mean_T1_D1_21(:,1)) transpose(mean_T1_D2_21(:,1)) transpose(mean_T1_D3_21(:,1)) subject_code(ind) ind transpose(mean_T1_D1_21(1,6)) transpose(mean_T1_D2_21(1,6)) transpose(mean_T1_D3_21(1,6)) h r];
                 T(ind,:) = [transpose(mean_T1_D1_21(:,1)) transpose(mean_T1_D2_21(:,1)) transpose(mean_T1_D3_21(:,1)) subject_code(ind) ind transpose(mean_T1_D1_21(1,6)) transpose(mean_T1_D2_21(1,6)) transpose(mean_T1_D3_21(1,6)) h r];
                     
                 % specify datasets that contain NaN values in QSM of layers #4 - #19:
                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                 tmp = [mean_T1_D1_21(:,1) mean_T1_D2_21(:,1) mean_T1_D3_21(:,1)];
                      
                 subjects_young(ind_y) = subject(ind);
                 color_y = sscanf(colorspec{ind},'%2x%2x%2x',[1 3])/255;

                 if (ind==1 & strcmp(hemi{h},'right')) || (ind == 2 & strcmp(hemi{h},'left'))
                    color = 'r'
                 else
                     color = 'k'
                 end

                 figure(ind); % individual color
                 for f=1:n_fingers
                    subplot(1,7,f); %f
                    plot(transpose(tmp(1:1:end,f)),'LineWidth',2.0,'Color',color);
                    %plot(transpose(tmp(1:1:end,f)),'LineWidth',2.0);
                    hold on;

                    %subplot(1,7,f+3); %f+3
                    %plot(abs(gradient(transpose(tmp(1:1:end,f)))),'LineWidth',2.0,'Color',color);
                    %plot(abs(gradient(transpose(tmp(1:1:end,f)))),'LineWidth',2.0);
                    %hold on;
                 end
        
                 cd (result_dir);
  
                 ind_y = ind_y + 1;
                    
             end

        end

    end

    saveas(figure(ind),sprintf('%s_pQSM_mFaceHandFoot_rotated_%s_70_minROI.eps',subject{ind},res{r}),'epsc');
    saveas(figure(ind),sprintf('%s_pQSM_mFaceHandFoot_rotated_%s_70_minROI.fig',subject{ind},res{r}));
    saveas(figure(ind),sprintf('%s_pQSM_mFaceHandFoot_rotated_%s_70_minROI.jpeg',subject{ind},res{r}));
            
end

cd (result_dir);

%save('allsubjects_qT1_mFaceHandFoot_short_resampled.mat','table_allsubjects_fingers_short');
%save('allsubjects_qT1_mFaceHandFoot_long_resampled.mat','table_allsubjects_fingers_long');

%csvwrite('allsubjects_qT1_mFaceHandFoot_short_resampled.csv',table_allsubjects_fingers_short);
%csvwrite('allsubjects_qT1_mFaceHandFoot_long_resampled.csv',table_allsubjects_fingers_long);
