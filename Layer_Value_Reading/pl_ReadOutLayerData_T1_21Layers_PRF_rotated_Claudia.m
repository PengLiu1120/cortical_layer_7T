%% Read out layer data
% this script reads in the finger maps, binarizes the finger maps, and multiplies the finger maps with the T1-map, and reads out data

%% ........................................................................Tidy up
clear all
close all
clc

%% ........................................................................Set defaults

% .........................................................................Subject definitions
subject = {'frj712' 'gxo876' 'hby152' 'ijt563' 'kdy341' 'lpr469' 'nhm378' 'oms448' 'qet940' 'qxo538' 'unk742'};
% subject = {'ajz367' 'bkn792' 'bmg520' 'cxc075' 'czg996' 'ggp057' 'gph998' 'iwq192' 'llh150' 'sst050'};

% .........................................................................File definitions
T1file_nameSurface = '*_surf_2-2_inf__all_layers.vtk';
T1file_nameSurfaceFolder = {'*AAAAAA'};

function_nameSurface = {'*_d1_smoothdata.vtk'};
function_nameSurfaceFolder = {'*BAAA'};

exp = {'exp-0005' 'exp-0008' 'exp-0009' 'exp-0011' 'exp-0013' 'exp-0014' 'exp-0015' 'exp-0016' 'exp-0017' 'exp-0018' 'exp-0020'};
% exp = {'exp-0000' 'exp-0001' 'exp-0002' 'exp-0003' 'exp-0004' 'exp-0006' 'exp-0007' 'exp-0010' 'exp-0012' 'exp-0019'};

% .........................................................................Path definitions
% .........................................................................RootDir
Rootdir = '/Volumes/LayerPRF/LayerMapping';
% .........................................................................T1 file folder
datadir = '04_LayerMapping'; 
% .........................................................................Functional file folder
datadir_func = '05_LayerDefining'; 

result_dir = '/Volumes/LayerPRF/LayerMapping/06_Derivative';

% .........................................................................Specifiy number of finger maps
n_fingers = 1;

table_allsubjects_hand_short = [];

n_layers = 21;
start_layer = 1;
end_layer = 21;

%% ........................................................................Start loop
% .........................................................................extract T1 values and finger maps for all subjects (one iteration per subject subject)
% .........................................................................store and plot data
for ind = 1:length(exp)

    Currexp = exp{ind};
    
    clear('T1_layers_area3b')
    clear('func_tval_finger')
    clear('finger_winner')
    
    funcfile_folder = fullfile(Rootdir, datadir_func, Currexp);
    T1file_folder = fullfile(Rootdir, datadir, Currexp);
    
    % .....................................................................Extact T1 values at different cortical depths
    % .....................................................................Go to T1file folder
    subjectdir = T1file_folder;
    cd (subjectdir)
    
    T1fileinfo = dir(T1file_nameSurfaceFolder{1});
    T1file = T1fileinfo(1).name;
    cd (T1file)
    cd ('SurfaceMeshMapping')
    T1fileinfo =  dir(T1file_nameSurface);
    T1file = T1fileinfo(1).name;
    
    [vertex,face,T1values,header1,header2,header3] = read_vtk_all(T1file);
    
    % .....................................................................Store T1 values
    T1_layers_area3b = T1values;
    
    % .....................................................................Start loop to extract prfCL values from functional maps:
    % .....................................................................read out prfCL values for all fingers (one iteration per finger)
    for ind3 = 1:n_fingers
        
        % .................................................................Go to data folder
        subjectdir = funcfile_folder;
        cd (subjectdir);
        
        funcfolderinfo =  dir(function_nameSurfaceFolder{ind3});
        funcfolder = funcfolderinfo(1).name;
        cd (funcfolder);
        cd ('SmoothSurfaceMeshData');
        funcfileinfo =  dir(function_nameSurface{ind3});
        funcfile = funcfileinfo(1).name;
        
        % .................................................................Read in surface file in vtk format
        [vertex,face,tval,header1,header2,header3] = read_vtk(funcfile);
        
        % .................................................................Store t-values to variable
        func_tval_finger(ind3,1:length(tval)) = tval;
        
    end
    
    % .....................................................................Threshold data at 0.1 to exclude smaller values by setting them to zero
    func_tval_finger(func_tval_finger<0.1) = 0;
    
    % .....................................................................Binarize finger maps & find overlaps by finding sum of columns > 1
    % .....................................................................(reflecting that vertex is occupied by more than one finger
    func_tval_finger_bin = func_tval_finger;
    func_tval_finger_bin(func_tval_finger_bin>0) = 1;
    
    cd (result_dir);
    
    % .....................................................................Multiply binarized ROI mask with T1 values
    T1_ROI = T1_layers_area3b.*func_tval_finger_bin(1,:);
    
    % .....................................................................Write data to matfile
    save(sprintf('%s_T1_ROI.mat',subject{ind}),'T1_ROI');
    
    % .....................................................................Replace zero with NaN
    T1_ROI(T1_ROI==0) = NaN;
    
    % .....................................................................Calculate mean across vertices per layer (one value per layer)
    [rownum,colnum]=size(T1_ROI);
    
    % .....................................................................Iterate over all layers
    for ind5 = 1:rownum
        
        % .................................................................Generates column vector with 4 rows (one row per layer containing average T1 value in the first column
        mean_T1_D1(ind5,1) = nanmean(T1_ROI(ind5,1:end));
        % .................................................................Generates columns vector with 4 rows containing finger number
        mean_T1_D1(ind5,2) = 1; 
        % .................................................................Generates column vector with 4 rows containing layer number
        mean_T1_D1(ind5,3) = ind5; 
        mean_T1_D1(ind5,4) = ind;
    end
    
    % .....................................................................Write data to matfile
    save(sprintf('%s_T1_ROI_mean.mat',subject{ind}),'mean_T1_D1');
    
    % .....................................................................Add all data to sample table
    table_allsubjects_hand_short (ind,:) = [transpose(mean_T1_D1(:,1)) ind];
    
end

cd (result_dir);

% .........................................................................Save figures & data
save('allsubjects_T1_young_hand_short.mat','table_allsubjects_hand_short');
% save('allsubjects_T1_old_hand_short.mat','table_allsubjects_hand_short');