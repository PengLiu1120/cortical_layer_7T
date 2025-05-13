%% Read out Fourier layer data
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
RootDir = '/media/pliu/LayerPRF/LayerMapping';

% .........................................................................Layer file folder
DataDir = '07_FunctionalMapping/D3_i_o';

% .........................................................................Results folder
ResultDir = '08_LayerExtraction/Fourier/D3';

%% ........................................................................Set defaults
% .........................................................................Specify subjects
Subjects = {'ajz367' 'bkn792' 'bmg520' 'cxc075' 'czg996' 'ggp057' 'gph998' 'iwq192' 'llh150' 'sst050'};
% 'frj712' 'gxo876' 'hby152' 'ijt563' 'kdy341' 'lpr469' 'nhm378' 'oms448' 'qet940' 'qxo538' 'unk742'
% 'ajz367' 'bkn792' 'bmg520' 'cxc075' 'czg996' 'ggp057' 'gph998' 'iwq192' 'llh150' 'sst050'

% .........................................................................Specify condtions, D2, D3 and D2+D3
Conditions = {'D2' 'D3' 'D2+D3'};

% .........................................................................Specify superficial layer results
Superficial = '*_surf_2-2_inf__S1_16-20_smoothdata.vtk';
Superficial_folder = '*AAAACAAA';

% .........................................................................Specify middle layer results
Middle = '*_surf_2-2_inf__S1_10-15_smoothdata.vtk';
Middle_folder = '*AAAABAAA';

% .........................................................................Specify deep layer results
Deep = '*__surf_2-2_inf__S1_2-9_smoothdata.vtk';
Deep_folder = '*AAAAAAAA';

% .........................................................................Specify exp
exp = {'exp-0000' 'exp-0001' 'exp-0002' 'exp-0003' 'exp-0004' 'exp-0006' 'exp-0007' 'exp-0010' 'exp-0012' 'exp-0019'};
% 'exp-0005' 'exp-0008' 'exp-0009' 'exp-0011' 'exp-0013' 'exp-0014' 'exp-0015' 'exp-0016' 'exp-0017' 'exp-0018' 'exp-0020'
% 'exp-0000' 'exp-0001' 'exp-0002' 'exp-0003' 'exp-0004' 'exp-0006' 'exp-0007' 'exp-0010' 'exp-0012' 'exp-0019'

%% ........................................................................Extract values for superficial layer from functional mapping results
% .........................................................................Layer 17 - 21
Superficial_Num = 5;

% .........................................................................Start loop
for i_sub=1:size(Subjects, 2)
    
    CurrSubj = Subjects{i_sub};

    Currexp = exp{i_sub};
    
    DataPath = fullfile(RootDir, DataDir, Currexp);
    cd(DataPath);
    
    % .................................................................Read superficial layer values
    Superficialfolder_info = dir(Superficial_folder);
    Superficial_file_folder = Superficialfolder_info(1).name;
    cd (Superficial_file_folder)
    cd ('SmoothSurfaceMeshData')
    Superficial_info =  dir(Superficial);
    Superficial_file = Superficial_info(1).name;
    
    [vertex,face,tval_superficial,header1,header2,header3] = read_vtk(Superficial_file);
    
    % .................................................................Store readout values from superficial layer
    superficial_layers = tval_superficial;
    
    ResultPath = fullfile(RootDir, ResultDir, CurrSubj);
    mkdir(ResultPath);
    cd(ResultPath);
    
    Superficial_Result = [CurrSubj '_Fourier_i_Superficial'];
    save(Superficial_Result, 'superficial_layers');
    
    % .................................................................Replace zero with NaN
    superficial_layers(superficial_layers==0) = NaN;
    superficial_layers(superficial_layers<=0.0001) = NaN;
    
    % .................................................................Calculate mean across vertices per layer (one value per layer)
    [rownum,colnum]=size(superficial_layers);
    
    % .................................................................Iterate over all layers
    for i_row = 1:rownum
        
        % .............................................................Generates column vector with 4 rows (one row per layer containing average T1 value in the first column
        mean_superficial(i_row,1) = nanmean(superficial_layers(i_row,1:end));
        % .............................................................Generates column vector with 4 rows containing layer number
        mean_superficial(i_row,2) = i_row;
        mean_superficial(i_row,3) = i_sub;
        
    end
    
    Superficial_Result_All = [CurrSubj '_Fourier_i_Superficial_layers'];
    save(Superficial_Result_All, 'mean_superficial');
    
end

%% ........................................................................Extract values for middle layer from functional mapping results
% .........................................................................Layer 11 - 16
Middle_Num = 6;

% .........................................................................Start loop
for i_sub=1:size(Subjects, 2)
    
    CurrSubj = Subjects{i_sub};
    
    Currexp = exp{i_sub};
    
    DataPath = fullfile(RootDir, DataDir, Currexp);
    cd(DataPath);
    
    % .....................................................................Read superficial layer values
    Middlefolder_info = dir(Middle_folder);
    Middle_file_folder = Middlefolder_info(1).name;
    cd (Middle_file_folder)
    cd ('SmoothSurfaceMeshData')
    Middle_info =  dir(Middle);
    Middle_file = Middle_info(1).name;
    
    [vertex,face,tval_middle,header1,header2,header3] = read_vtk(Middle_file);
    
    % .....................................................................Store readout values from superficial layer
    middle_layers = tval_middle;
    
    ResultPath = fullfile(RootDir, ResultDir, CurrSubj);
    cd(ResultPath);
    
    Middle_Result = [CurrSubj '_Fourier_i_Middle'];
    save(Middle_Result, 'middle_layers');
    
    % .................................................................Replace zero with NaN
    middle_layers(middle_layers==0) = NaN;
    middle_layers(middle_layers<=0.0001) = NaN;
    
    % .................................................................Calculate mean across vertices per layer (one value per layer)
    [rownum,colnum]=size(middle_layers);
    
    % .................................................................Iterate over all layers
    for i_row = 1:rownum
        
        % .............................................................Generates column vector with 4 rows (one row per layer containing average T1 value in the first column
        mean_middle(i_row,1) = nanmean(middle_layers(i_row,1:end));
        % .............................................................Generates column vector with 4 rows containing layer number
        mean_middle(i_row,2) = i_row;
        mean_middle(i_row,3) = i_sub;
        
    end
    
    Middle_Result_All = [CurrSubj '_Fourier_i_Middle_layers'];
    save(Middle_Result_All, 'mean_middle');
    
end

%% ........................................................................Extract values for deep layer from functional mapping results
% .........................................................................Layer 3 - 10
Deep_Num = 8;

% .........................................................................Start loop
for i_sub=1:size(Subjects, 2)
    
    CurrSubj = Subjects{i_sub};
    
    Currexp = exp{i_sub};
    
    DataPath = fullfile(RootDir, DataDir, Currexp);
    cd(DataPath);
    
    % .....................................................................Read deep layer values
    Deepfolder_info = dir(Deep_folder);
    Deep_file_folder = Deepfolder_info(1).name;
    cd (Deep_file_folder)
    cd ('SmoothSurfaceMeshData')
    Deep_info =  dir(Deep);
    Deep_file = Deep_info(1).name;
    
    [vertex,face,tval_deep,header1,header2,header3] = read_vtk(Deep_file);
    
    % .....................................................................Store readout values from deep layer
    deep_layers = tval_deep;
    
    ResultPath = fullfile(RootDir, ResultDir, CurrSubj);
    cd(ResultPath);
    
    Deep_Result = [CurrSubj '_Fourier_i_Deep'];
    save(Deep_Result, 'deep_layers');
    
    % .................................................................Replace zero with NaN
    deep_layers(deep_layers==0) = NaN;
    deep_layers(deep_layers<=0.0001) = NaN;
    
    % .................................................................Calculate mean across vertices per layer (one value per layer)
    [rownum,colnum]=size(deep_layers);
    
    % .................................................................Iterate over all layers
    for i_row = 1:rownum
        
        % .............................................................Generates column vector with 4 rows (one row per layer containing average T1 value in the first column
        mean_deep(i_row,1) = nanmean(deep_layers(i_row,1:end));
        % .............................................................Generates column vector with 4 rows containing layer number
        mean_deep(i_row,2) = i_row;
        mean_deep(i_row,3) = i_sub;
        
    end
    
    Deep_Result_All = [CurrSubj '_Fourier_i_Deep_layers'];
    save(Deep_Result_All, 'mean_deep');
    
end