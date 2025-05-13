%% Read out layer-specific resting state data
% .........................................................................
% This script reads out resting state data at layer-specific for 19 layers.
% .........................................................................
% Written by P.Liu
% Optimized by P.Liu
% Email: peng.liu@med.ovgu.de
% Last modified by P.Liu 13 May 2023
%% ........................................................................Tidy up
clear all
close all
clc

%% ........................................................................Set paths
% .........................................................................Specify RootDir
RootDir = '/Volumes/IKND/AG_Kuehn/Peng/LayerPRF/LayerMapping';

% .........................................................................Layer file folder
DataDir = '07_FunctionalMapping/Resting_State_mean';

% .........................................................................Results folder
ResultDir = '08_LayerExtraction/Resting_State_mean';

%% ........................................................................Set defaults
% .........................................................................Specify subjects
Subjects = {'frj712' 'gxo876' 'hby152' 'ijt563' 'kdy341' 'lpr469' 'nhm378' 'oms448' 'qet940' 'qxo538' 'unk742' 'ajz367' 'bkn792' 'bmg520' 'cxc075' 'czg996' 'ggp057' 'gph998' 'iwq192' 'llh150' 'sst050'};

% .........................................................................Specify exp
exp = {'exp-0005' 'exp-0008' 'exp-0009' 'exp-0011' 'exp-0013' 'exp-0014' 'exp-0015' 'exp-0016' 'exp-0017' 'exp-0018' 'exp-0020' 'exp-0000' 'exp-0001' 'exp-0002' 'exp-0003' 'exp-0004' 'exp-0006' 'exp-0007' 'exp-0010' 'exp-0012' 'exp-0019'};

% .........................................................................Specify layer folders
folders = {'A' 'B' 'C' 'D' 'E' 'F' 'G' 'H' 'I' 'J' 'K' 'L' 'M' 'N' 'O' 'P' 'Q' 'R' 'S'};

% .........................................................................Specify layer numbers
numbers = {'2' '3' '4' '5' '6' '7' '8' '9' '10' '11' '12' '13' '14' '15' '16' '17' '18' '19' '20'};

%% ........................................................................Extract values for each layer from resting state mapping results
for i_sub=1:size(Subjects, 2)
    
    CurrSubj = Subjects{i_sub};
    
    Currexp = exp{i_sub};
    
    DataPath = fullfile(RootDir, DataDir, Currexp);
    cd(DataPath);
    
    tval_layers = [];
    mean_layers = [];
    
    for num=1:19
        
        CurrLayer = numbers{num};
        CurrFolder = folders{num};
        
        layer = ['*_surf_2-2_inf__S1_' CurrLayer '_smoothdata.vtk'];
        layer_folder = ['*AAAA' CurrFolder 'AA'];
        
        layer_folder_info = dir(layer_folder);
        layer_file_folder = layer_folder_info(1).name;
        
        cd(layer_file_folder);
        cd('SmoothSurfaceMeshData');
        
        layer_info = dir(layer);
        layer_file = layer_info(1).name;
        
        [vertex,face,tval_layer,header1,header2,header3] = read_vtk(layer_file);
        
        % .................................................................Store readout values from each layer
        layers_3b = tval_layer;
        tval_layers = [tval_layers tval_layer];
        
        % .................................................................Replace zero with NaN
        layers_3b(layers_3b==0) = NaN;
        layers_3b(layers_3b<=0.0001) = NaN;
        
        mean_layer = nanmean(layers_3b);
        mean_layers = [mean_layers; mean_layer];
        
        cd(DataPath);
        
    end
    
    ResultPath = fullfile(RootDir, ResultDir, CurrSubj);
    cd(ResultPath);
    
    layer_result = [CurrSubj '_resting_state_layer_3b'];
    save(layer_result, 'tval_layers');
    
    All_Layers_Result = [CurrSubj '_resting_state_all_layers_3b'];
    save(All_Layers_Result, 'mean_layers');
    
end