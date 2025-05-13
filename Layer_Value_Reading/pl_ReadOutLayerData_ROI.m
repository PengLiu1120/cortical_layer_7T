%% Read out layer data
% this script reads in the finger maps, binarizes the finger maps and
% multiplies the finger maps with the ROI mask, and reads out data
% .........................................................................
% Superficial layer
% Deep layer
% .........................................................................
% Written by P.Liu
% Optimized by P.Liu
% Email: peng.liu@med.ovgu.de
% Last modified by P.Liu 28 Oct 2022
%% ........................................................................Tidy up
clear all
close all
clc

%% ........................................................................Set paths
% .........................................................................Specify RootDir
RootDir = '/Users/pliu/Documents/LayerPRF/LayerMapping';

% .........................................................................ROI folder
ROIDir = '05_LayerDefining'; 

% .........................................................................Results folder
ResultDir = '08_Results';

%% ........................................................................Set defaults
% .........................................................................Specify subjects
Subjects = {'frj712' 'lpr469' 'qet940' 'qxo538' 'oms448' 'hby152' 'ijt563'};
% 'frj712' 'lpr469' 'qet940' 'qxo538' 'oms448' 'hby152' 'ijt563'

% .........................................................................Specify ROI mask
ROI = '*_d1_smoothdata.vtk';
ROI_folder = '*BAAA';

% .........................................................................Specify exp
ROI_exp = 'exp-0000';

%% ........................................................................Extract ROI mask
% .........................................................................Start loop
for i_sub=1:size(Subjects, 2)
    
    CurrSubj = Subjects{i_sub};
    
    % .....................................................................Read ROI mask
    ROIPath = fullfile(RootDir, ROIDir, CurrSubj, ROI_exp);
    cd(ROIPath);
    
    ROIfolder_info =  dir(ROI_folder);
    ROI_file_folder = ROIfolder_info(1).name;
    cd (ROI_file_folder);
    cd ('SmoothSurfaceMeshData');
    ROI_info =  dir(ROI);
    ROI_file = ROI_info(1).name;
    
    % .....................................................................Read in surface file in vtk format
    [vertex,face,tval,header1,header2,header3] = read_vtk(ROI_file);
    
    % .....................................................................Store t-values to variable
    ROI_tval_finger(1,1:length(tval)) = tval;
    
    
    % .....................................................................Threshold data at 0.1 to exclude smaller values by setting them to zero
    ROI_tval_finger(ROI_tval_finger<0.1) = 0;
    
    % .....................................................................Binarize finger maps & find overlaps by finding sum of columns > 1
    % .....................................................................(reflecting that vertex is occupied by more than one finger
    ROI_tval_finger_bin = ROI_tval_finger;
    ROI_tval_finger_bin(ROI_tval_finger_bin>0) = 1;
    
    ResultPath = fullfile(RootDir, ResultDir);
    cd(ResultPath);
    
    ROI_FileName = ['ROI_' CurrSubj];
    save(ROI_FileName, 'ROI_tval_finger_bin');
    
end