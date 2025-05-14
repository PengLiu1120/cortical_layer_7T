% this script reads in the finger maps, binarizes the finger maps, and multiplies the
% finger maps with the T1-map, and reads out data

clear;

%% General Specifications
%%%%%%%%%%%%%%%%%%%%%%%%%
fileID = fopen('/home/juli/projects/onehander/data/paths/all_IDs.txt','r');
formatSpec = '%3s';

A = textscan(fileID,formatSpec)

subject = A{1,1}'

%subject_folder = {'00','01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','20','21','22','23','24','25','26','27','28','29','30','31','32','33','34','35','36','37','38','39','40'};

subject_code = [3 3]; % 1=young, 2=old, 3=onehander

% for header information only
% S1 data
T1file_nameSurface = '*_surf_2-2_inf__S1_3-6.vtk';
headerdir = '/home/juli/projects/onehander/data';

function_nameSurface = {'*_hand_S1_smoothdata.vtk','*_face_S1_smoothdata.vtk','*_foot_S1_smoothdata.vtk'};

% definitions
% for laptop
datadir = '/home/juli/projects/onehander/data';
datadir_func = '/home/juli/projects/onehander/data';

script_dir = '/home/juli/projects/onehander/scripts';
result_dir = '/home/juli/projects/onehander/results/ReadOutLayerData_T1_HandFaceFoot';
out_dir = '/home/juli/projects/onehander/results/ReadOutLayerData_QSM_HandFaceFoot';

layer_def = '/home/juli/projects/onehander/results/extract_Anatomical_layers';

addpath(script_dir)

% specifiy number of functional maps
n_fingers = 3;

n_layers = 3;

%file_ext_loc = {'hand','face','foot'};

file_ext_loc = {'inner','middle','outer','inner','middle','outer','inner','middle','outer','inner','middle','outer'};

qsm_type = {'aQSM','aQSM','aQSM','pQSM','pQSM','pQSM','nQSM','nQSM','nQSM','sQSM','sQSM','sQSM'};

hemi = {'left','right'};

%res = {'orig','resampled'};
res = {'orig'};

table_allsubjects_fingers_long_pos = [];
table_allsubjects_fingers_short_pos = [];

table_allsubjects_fingers_long_21 = [];
table_allsubjects_fingers_short_21 = [];

table_allsubjects_fingers_long_neg = [];
table_allsubjects_fingers_short_neg = [];

table_allsubjects_fingers_long_21_neg = [];
table_allsubjects_fingers_short_21_neg = [];

table_allsubjects_fingers_long_abs = [];
table_allsubjects_fingers_short_abs = [];

table_allsubjects_fingers_long_21_abs = [];
table_allsubjects_fingers_short_21_abs = [];

%% start loop to: extract T1 values and t-values for all subjects (one iteration per subject subject), store and plot data
for ind = 2:length(subject)

    for h = 1:length(hemi)  

        if ind == 2 && h == 1
            function_nameSurfaceFolder = {'*ECAA','*EDAA','*EHAA'};
        else
            function_nameSurfaceFolder = {'*ECAA','*EDAA','*EEAA'};
        end
    
        for r = 1:length(res)
    
         clear('T1_layers_area3b')
         clear('T1_layers_area3b_bin')
         clear('func_tval_finger')
         clear('finger_winner')
         clear('finger_winner_S1')

              % load layer definition
         if h == 1
             load(sprintf('%s/%s/idx_innerLayer_l_a_70_minROI.mat',layer_def,subject{ind}));
             load(sprintf('%s/%s/idx_middleLayer_l_a_70_minROI.mat',layer_def,subject{ind}));
             load(sprintf('%s/%s/idx_outerLayer_l_a_70_minROI.mat',layer_def,subject{ind}));
             ind_inner = ind_inner;
             ind_middle = ind_middle;
             ind_outer = ind_outer;
         elseif h == 2
             load(sprintf('%s/%s/idx_innerLayer_r_na_70_minROI.mat',layer_def,subject{ind}));
             load(sprintf('%s/%s/idx_middleLayer_r_na_70_minROI.mat',layer_def,subject{ind}));
             load(sprintf('%s/%s/idx_outerLayer_r_na_70_minROI.mat',layer_def,subject{ind}));
             ind_inner = ind_inner; %ind_inner_o
             ind_middle = ind_middle; %ind_middle_o
             ind_outer = ind_outer; %ind_outer_o
         end

        inner=(diff(ind_inner(1,1:2))+1)/2
        middle=(diff(ind_middle(1,1:2))+1)/2
        outer=(diff(ind_outer(1,1:2))+1)/2
    
        funcfile_folder = 'exp-0000'

        % go to T1file folder and load extracted qT1 values for 21 layers
        cd (sprintf('%s/%s/processed/output/ReadOutLayerData_QSM_21Layers_withfoot',datadir, subject{ind}));
        load(sprintf('%s_pQSM_BA3b_21layers_%s.mat',subject{ind},hemi{h}));
        load(sprintf('%s_nQSM_BA3b_21layers_%s.mat',subject{ind},hemi{h}));
        load(sprintf('%s_aQSM_BA3b_21layers_%s.mat',subject{ind},hemi{h}));
        load(sprintf('%s_sQSM_BA3b_21layers_%s.mat',subject{ind},hemi{h}));
     
        area3b = T1_layers_area3b;
        area3b_neg = T1_layers_area3b_neg;
        area3b_abs = T1_layers_area3b_abs;
        sig_qsm_all_layers_area3b = T1_area3b;
        
        T1_layers_area3b = [];
    
        % replace zero with NaN
        T1_layers_area3b_NaN = area3b;
        T1_layers_area3b_NaN(T1_layers_area3b_NaN==0) = NaN;

        T1_layers_area3b_NaN_neg = area3b_neg;
        T1_layers_area3b_NaN_neg(T1_layers_area3b_NaN_neg==0) = NaN;

        T1_layers_area3b_NaN_abs = area3b_abs;
        T1_layers_area3b_NaN_abs(T1_layers_area3b_NaN_abs==0) = NaN;
           
        % % generate binary S1 mask
        T1_layers_area3b_bin = T1_layers_area3b_abs(10,:);
        T1_layers_area3b_bin(isnan(T1_layers_area3b_bin)) = 0;
        T1_layers_area3b_bin(T1_layers_area3b_bin>0) = 1;

        qsm_layers_area3b = [];
        qsm_layers_area3b = zeros((3)*3,length(area3b_abs));
                  
        % average qT1 values into anatomically relevant layer compartments
        qsm_layers_area3b(1,find(sum(area3b_abs(ind_inner(1,1):ind_inner(1,2),:)>0,1)>inner)) = nanmean(area3b_abs(ind_inner(1,1):ind_inner(1,2),find(sum(area3b_abs(ind_inner(1,1):ind_inner(1,2),:)>0,1)>inner)),1);
        qsm_layers_area3b(2,find(sum(area3b_abs(ind_middle(1,1):ind_middle(1,2),:)>0,1)>middle)) = nanmean(area3b_abs(ind_middle(1,1):ind_middle(1,2),find(sum(area3b_abs(ind_middle(1,1):ind_middle(1,2),:)>0,1)>middle)),1);
        qsm_layers_area3b(3,find(sum(area3b_abs(ind_outer(1,1):ind_outer(1,2),:)>0,1)>outer)) = nanmean(area3b_abs(ind_outer(1,1):ind_outer(1,2),find(sum(area3b_abs(ind_outer(1,1):ind_outer(1,2),:)>0,1)>outer)),1);
        %qsm_layers_area3b(4,find(sum(area3b_abs(15:18,:)>0,1)>2)) = nanmean(area3b_abs(15:18,find(sum(area3b_abs(15:18,:)>0,1)>2)),1);
        qsm_layers_area3b(4,find(sum(area3b(ind_inner(1,1):ind_inner(1,2),:)>0,1)>inner)) = nanmean(area3b(ind_inner(1,1):ind_inner(1,2),find(sum(area3b(ind_inner(1,1):ind_inner(1,2),:)>0,1)>inner)),1);
        qsm_layers_area3b(5,find(sum(area3b(ind_middle(1,1):ind_middle(1,2),:)>0,1)>middle)) = nanmean(area3b(ind_middle(1,1):ind_middle(1,2),find(sum(area3b(ind_middle(1,1):ind_middle(1,2),:)>0,1)>middle)),1);
        qsm_layers_area3b(6,find(sum(area3b(ind_outer(1,1):ind_outer(1,2),:)>0,1)>outer)) = nanmean(area3b(ind_outer(1,1):ind_outer(1,2),find(sum(area3b(ind_outer(1,1):ind_outer(1,2),:)>0,1)>outer)),1);
        %qsm_layers_area3b(8,find(sum(area3b(15:18,:)>0,1)>2)) = nanmean(area3b(15:18,find(sum(area3b(15:18,:)>0,1)>2)),1);
        qsm_layers_area3b(7,find(sum(area3b_neg(ind_inner(1,1):ind_inner(1,2),:)<0,1)>inner)) = nanmean(area3b_neg(ind_inner(1,1):ind_inner(1,2),find(sum(area3b_neg(ind_inner(1,1):ind_inner(1,2),:)<0,1)>inner)),1);
        qsm_layers_area3b(8,find(sum(area3b_neg(ind_middle(1,1):ind_middle(1,2),:)<0,1)>middle)) = nanmean(area3b_neg(ind_middle(1,1):ind_middle(1,2),find(sum(area3b_neg(ind_middle(1,1):ind_middle(1,2),:)<0,1)>middle)),1);
        qsm_layers_area3b(9,find(sum(area3b_neg(ind_outer(1,1):ind_outer(1,2),:)<0,1)>outer)) = nanmean(area3b_neg(ind_outer(1,1):ind_outer(1,2),find(sum(area3b_neg(ind_outer(1,1):ind_outer(1,2),:)<0,1)>outer)),1);
        %qsm_layers_area3b(12,find(sum(area3b_neg(15:18,:)<0,1)>2)) = nanmean(area3b_neg(15:18,find(sum(area3b_neg(15:18,:)<0,1)>2)),1);
        qsm_layers_area3b(10,find(sum(sig_qsm_all_layers_area3b(ind_inner(1,1):ind_inner(1,2),:)~=0,1)>inner)) = nanmean(sig_qsm_all_layers_area3b(ind_inner(1,1):ind_inner(1,2),find(sum(sig_qsm_all_layers_area3b(ind_inner(1,1):ind_inner(1,2),:)~=0,1)>inner)),1);
        qsm_layers_area3b(11,find(sum(sig_qsm_all_layers_area3b(ind_middle(1,1):ind_middle(1,2),:)~=0,1)>middle)) = nanmean(sig_qsm_all_layers_area3b(ind_middle(1,1):ind_middle(1,2),find(sum(sig_qsm_all_layers_area3b(ind_middle(1,1):ind_middle(1,2),:)~=0,1)>middle)),1);
        qsm_layers_area3b(12,find(sum(sig_qsm_all_layers_area3b(ind_outer(1,1):ind_outer(1,2),:)~=0,1)>outer)) = nanmean(sig_qsm_all_layers_area3b(ind_outer(1,1):ind_outer(1,2),find(sum(sig_qsm_all_layers_area3b(ind_outer(1,1):ind_outer(1,2),:)~=0,1)>outer)),1);

%         % replace zero with NaN
%         qsm_layers_area3b_NaN = qsm_layers_area3b;
%         qsm_layers_area3b_NaN(qsm_layers_area3b_NaN==0) = NaN;
% 
%         % Calculate mean value per layer compartment     
%         for ind2 = 1:n_layers
%            % Calculate mean      
%            mean_QSM_layers_3b(ind,ind2) = nanmean(qsm_layers_area3b_NaN(ind2,1:end)); %absolute values, reflecting level of mineralisation
%            mean_pos_QSM_layers_3b(ind,ind2) = nanmean(qsm_layers_area3b_NaN(ind2+n_layers,1:end)); %positive values only, reflection iron content
%            mean_neg_QSM_layers_3b(ind,ind2) = nanmean(qsm_layers_area3b_NaN(ind2+(n_layers*2),1:end)); %negative values only, refelecting myelin/calcium content 
%            %mean_sig_QSM_layers_3b(qsm_ind,ind2) = nanmean(qsm_layers_area3b_NaN(ind2+(n_layers*3),1:end)); %negative values only, refelecting myelin/calcium content 
%         end
% 
%         % write data to matfile
%         cd(result_dir);
%         save(sprintf('%s_qsm_Layers_70_minROI.mat',subject{ind}),'mean_QSM_layers_3b','mean_pos_QSM_layers_3b','mean_neg_QSM_layers_3b'); %change file name to _Layers.mat for group-specific outcomes
% 
%         % write anatomically relevant layer maps to vtk files
%         qsm_layers_area3b(isnan(qsm_layers_area3b)) = 0;
% 
%         % start loop to extract values from functional maps: read out prfCL values for all fingers (one iteration per finger)
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for ind3 = 1%n_fingers

            % go to data folder
            cd (sprintf('%s/%s/processed/output/Mapping_hand_face_foot_FINAL_%s/%s',datadir_func,subject{ind},hemi{h},funcfile_folder));

            try
                funcfolderinfo =  dir(function_nameSurfaceFolder{ind3});
                funcfolder = funcfolderinfo(1).name;
                face_hand = 1;
            catch
                warning('No 1.5 mm localizer available. Trying next subject');
                face_hand = 0;
                ind3 = n_fingers+1;
            end

            if face_hand == 1
                cd (funcfolder);
                cd ('SmoothSurfaceMeshData');
                funcfileinfo = dir(function_nameSurface{ind3});
                funcfile = funcfileinfo(1).name

                % read in surface file in vtk format
                [vertex,face,tval,header1,header2,header3] = read_vtk(funcfile);

                % store t-values to variable
                func_tval_finger(ind3,1:length(tval)) = tval.*T1_layers_area3b_bin;
            end

        end

        
        cd(out_dir)
        for i = 10:12
            filename_new = sprintf('%s_%s_%s_%s_%s.vtk',subject{ind},file_ext_loc{i},qsm_type{i},hemi{h},res{r}); 
            write_vtk(filename_new,header1,header2,header3,vertex,face,qsm_layers_area3b(i,:))
        end
% 
%         if face_hand == 1
% 
%          % winner takes it all: takes care of overlapping finger areas (introduced by remapping divided prf digit maps and smoothing), therefore only keep the maximum finger value if there are multiple values per voxel (multiple values per column), set all other finger values in the same row to zero 
%          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          % M = vector including maximum values of each column; I = vector including row indeces of maximum values
%          % M also gives mask for all five fingers/hand
%          [M,I] = max(func_tval_finger);
%          M_bin = M;
%          M_bin(M_bin~=0) = 1;
% 
%          % to separate fingers into different rows define winner array of zeros
%          finger_winner = zeros(n_fingers,length(M));
% 
%          % start loop to fill matrix of zeros with maximum prf var values extacted from the different finger maps
%          for ind5 = 1:length(M)
%             finger_winner(I(ind5),ind5) = M(ind5);
%          end
% 
%          length_winner = [];
% 
%          for win = 1:size(finger_winner,1)
%              tmp = finger_winner(win,:);
%              tmp = tmp(tmp>0);
%              length_winner(win) = length(tmp);
%          end
% 
%          for win = 1:size(finger_winner,1)
%               tmp_sort = sort(finger_winner(win,:));
%              %min(length_winner)
%              top90 = round(length_winner(win)*0.7)
%              top = tmp_sort(end-1788)
%              %top = tmp_sort(end-top90)
% 
%              tmp_win = finger_winner(win,:);
%              tmp_win(tmp_win<top) = 0;
% 
%              tmp = tmp_win;
%              tmp = tmp(tmp>0);
%              length(tmp)
% 
%              finger_winner(win,:) = tmp_win;
% 
%              finger_winner_zvals = vertex(3,finger_winner(win,:)>0);
% 
%              %min(finger_winner_zvals)
% 
%              finger_winner_zvals_mean = mean(finger_winner_zvals);
%              finger_winner_zvals_sd = std(finger_winner_zvals);
% 
%              %test = 'finger_winner_tmp'
% 
%              finger_winner_tmp = finger_winner(win,:);
%              finger_winner_tmp(1,vertex(3,:)>finger_winner_zvals_mean+(2*finger_winner_zvals_sd)) = 0;
%              finger_winner_tmp(1,vertex(3,:)<finger_winner_zvals_mean-(2*finger_winner_zvals_sd)) = 0;
% 
%              finger_winner(win,:) = finger_winner_tmp;
% 
%              %if ind3 == 3
%                 finger_winner_yvals = vertex(2,finger_winner(win,:)>0);
%                 finger_winner_yvals_mean = mean(finger_winner_yvals);
%                 finger_winner_yvals_sd = std(finger_winner_yvals);
% 
%                 finger_winner_tmp(1,vertex(2,:)>finger_winner_yvals_mean+(2*finger_winner_yvals_sd)) = 0;
%                 finger_winner_tmp(1,vertex(2,:)<finger_winner_yvals_mean-(2*finger_winner_yvals_sd)) = 0;
% 
%                 finger_winner(win,:) = finger_winner_tmp;
% 
%                 finger_winner_xvals = vertex(1,finger_winner(win,:)>0);
%                 finger_winner_xvals_mean = mean(finger_winner_xvals);
%                 finger_winner_xvals_sd = std(finger_winner_xvals);
% 
%                 finger_winner_tmp(1,vertex(1,:)>finger_winner_xvals_mean+(2*finger_winner_xvals_sd)) = 0;
%                 finger_winner_tmp(1,vertex(1,:)<finger_winner_xvals_mean-(2*finger_winner_xvals_sd)) = 0;
% 
%                 finger_winner(win,:) = finger_winner_tmp;
%              %end
% 
%              %min(vertex(2,finger_winner(win,:)>0))
%          end
% 
%          cd (out_dir);
% 
%          %binarize
%          finger_winner_bin = finger_winner;
%          finger_winner_bin(finger_winner_bin~=0) = 1;
% 
%          sum(finger_winner_bin(1,:)>0)
%          sum(finger_winner_bin(2,:)>0)
%          sum(finger_winner_bin(3,:)>0)
% 
%          % apply S1 mask to final finger maps
%          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          finger_winner_S1 = zeros(n_fingers,length(finger_winner));
% 
%          for i = 1:n_fingers
%             finger_winner_S1(i,:) = finger_winner(i,:);
% 
%             %filename_new = sprintf('%s_%s_win_thrs_%s_%s.vtk',subject{ind},file_ext_loc{i},hemi{h},res{r}); 
%             %write_vtk(filename_new,header1,header2,header3,vertex,face,finger_winner_S1(i,:))
%          end
% 
%          % write data to matfile
%          %save(sprintf('%s_body_part_win_S1_thrs_%s_%s.mat',subject{ind},hemi{h},res{r}),'finger_winner_S1');
% 
%          %binarize
%          finger_winner_S1_bin = finger_winner_S1;
%          finger_winner_S1_bin(finger_winner_S1_bin~=0) = 1;
% 
%          % multiply binarized finger masks with T1 values
%          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%          T1_D1_21 = area3b.*finger_winner_bin(1,:); %d1 21 layers
%          T1_D2_21 = area3b.*finger_winner_bin(2,:); %d2 21 layers
%          T1_D3_21 = area3b.*finger_winner_bin(3,:); %d2 21 layers
%          T1_all_21 = area3b.*M_bin; %hand (all fingers) %all 21 layers  
% 
%          T1_D1_21_neg = area3b_neg.*finger_winner_bin(1,:); %d1 21 layers
%          T1_D2_21_neg = area3b_neg.*finger_winner_bin(2,:); %d2 21 layers
%          T1_D3_21_neg = area3b_neg.*finger_winner_bin(3,:); %d2 21 layers
%          T1_all_21_neg = area3b_neg.*M_bin; %hand (all fingers) %all 21 layers 
% 
%          T1_D1_21_abs = area3b_abs.*finger_winner_bin(1,:); %d1 21 layers
%          T1_D2_21_abs = area3b_abs.*finger_winner_bin(2,:); %d2 21 layers
%          T1_D3_21_abs = area3b_abs.*finger_winner_bin(3,:); %d2 21 layers
%          T1_all_21_abs = area3b_abs.*M_bin; %hand (all fingers) %all 21 layers 
% 
%          T1_layers_area3b = qsm_layers_area3b
%          T1_D1 = T1_layers_area3b.*finger_winner_bin(1,:); %d1 DD layers
%          T1_D2 = T1_layers_area3b.*finger_winner_bin(2,:); %d2 DD layers
%          T1_D3 = T1_layers_area3b.*finger_winner_bin(3,:); %d2 DD layers
%          T1_all = T1_layers_area3b.*M_bin; %hand (all fingers) %all DD layers
% 
%          % calculate # of survived vertices for each finger map
%          vertices_d1_S1 = sum(T1_D1_21(10,:)>0)
%          vertices_d2_S1 = sum(T1_D2_21(10,:)>0)
%          vertices_d3_S1 = sum(T1_D3_21(10,:)>0)
% 
%          vertices_d1_S1_neg = sum(T1_D1_21_neg(10,:)<0)
%          vertices_d2_S1_neg = sum(T1_D2_21_neg(10,:)<0)
%          vertices_d3_S1_neg = sum(T1_D3_21_neg(10,:)<0)
% 
%          vertices_d1_S1_abs = sum(T1_D1_21_abs(10,:)>0)
%          vertices_d2_S1_abs = sum(T1_D2_21_abs(10,:)>0)
%          vertices_d3_S1_abs = sum(T1_D3_21_abs(10,:)>0)
% 
%          % write data to matfile
% 
%          save(sprintf('%s_pQSM_mhand_21layer_thrs_%s_%s_70_minROI.mat',subject{ind},hemi{h},res{r}),'T1_D1_21');
%          save(sprintf('%s_pQSM_mface_21layer_thrs_%s_%s_70_minROI.mat',subject{ind},hemi{h},res{r}),'T1_D2_21');
%          save(sprintf('%s_pQSM_mfoot_21layer_thrs_%s_%s_70_minROI.mat',subject{ind},hemi{h},res{r}),'T1_D3_21');
%          save(sprintf('%s_pQSM_mface_mhand_mfoot_21layer_thrs_%s_%s_70_minROI.mat',subject{ind},hemi{h},res{r}),'T1_all_21');
% % 
%          save(sprintf('%s_nQSM_mhand_21layer_thrs_%s_%s_70_minROI.mat',subject{ind},hemi{h},res{r}),'T1_D1_21_neg');
%          save(sprintf('%s_nQSM_mface_21layer_thrs_%s_%s_70_minROI.mat',subject{ind},hemi{h},res{r}),'T1_D2_21_neg');
%          save(sprintf('%s_nQSM_mfoot_21layer_thrs_%s_%s_70_minROI.mat',subject{ind},hemi{h},res{r}),'T1_D3_21_neg');
%          save(sprintf('%s_nQSM_mface_mhand_mfoot_21layer_thrs_%s_%s_70_minROI.mat',subject{ind},hemi{h},res{r}),'T1_all_21_neg');
% % 
%          save(sprintf('%s_aQSM_mhand_21layer_thrs_%s_%s_70_minROI.mat',subject{ind},hemi{h},res{r}),'T1_D1_21_abs');
%          save(sprintf('%s_aQSM_mface_21layer_thrs_%s_%s_70_minROI.mat',subject{ind},hemi{h},res{r}),'T1_D2_21_abs');
%          save(sprintf('%s_aQSM_mfoot_21layer_thrs_%s_%s_70_minROI.mat',subject{ind},hemi{h},res{r}),'T1_D3_21_abs');
%          save(sprintf('%s_aQSM_mface_mhand_mfoot_21layer_thrs_%s_%s_70_minROI.mat',subject{ind},hemi{h},res{r}),'T1_all_21_abs');
% 
%          T1_D1_tmp = T1_D1(1:3,:);
%          T1_D2_tmp = T1_D2(1:3,:);
%          T1_D3_tmp = T1_D3(1:3,:);
%          T1_all_tmp = T1_all(1:3,:);
%          save(sprintf('%s_aQSM_mhand_DDlayer_thrs_%s_%s_70_minROI.mat',subject{ind},hemi{h},res{r}),'T1_D1_tmp');
%          save(sprintf('%s_aQSM_mface_DDlayer_thrs_%s_%s_70_minROI.mat',subject{ind},hemi{h},res{r}),'T1_D2_tmp');
%          save(sprintf('%s_aQSM_mfoot_DDlayer_thrs_%s_%s_70_minROI.mat',subject{ind},hemi{h},res{r}),'T1_D3_tmp');
%          save(sprintf('%s_aQSM_mface_mhand_mfoot_DDlayer_thrs_%s_%s_70_minROI.mat',subject{ind},hemi{h},res{r}),'T1_all_tmp');
% 
%          T1_D1_tmp = T1_D1(4:6,:);
%          T1_D2_tmp = T1_D2(4:6,:);
%          T1_D3_tmp = T1_D3(4:6,:);
%          T1_all_tmp = T1_all(4:6,:);
%          save(sprintf('%s_pQSM_mhand_DDlayer_thrs_%s_%s_70_minROI.mat',subject{ind},hemi{h},res{r}),'T1_D1_tmp');
%          save(sprintf('%s_pQSM_mface_DDlayer_thrs_%s_%s_70_minROI.mat',subject{ind},hemi{h},res{r}),'T1_D2_tmp');
%          save(sprintf('%s_pQSM_mfoot_DDlayer_thrs_%s_%s_70_minROI.mat',subject{ind},hemi{h},res{r}),'T1_D3_tmp');
%          save(sprintf('%s_pQSM_mface_mhand_mfoot_DDlayer_thrs_%s_%s_70_minROI.mat',subject{ind},hemi{h},res{r}),'T1_all_tmp');
% 
%          T1_D1_tmp = T1_D1(7:9,:);
%          T1_D2_tmp = T1_D2(7:9,:);
%          T1_D3_tmp = T1_D3(7:9,:);
%          T1_all_tmp = T1_all(7:9,:);
%          save(sprintf('%s_nQSM_mhand_DDlayer_thrs_%s_%s_70_minROI.mat',subject{ind},hemi{h},res{r}),'T1_D1_tmp');
%          save(sprintf('%s_nQSM_mface_DDlayer_thrs_%s_%s_70_minROI.mat',subject{ind},hemi{h},res{r}),'T1_D2_tmp');
%          save(sprintf('%s_nQSM_mfoot_DDlayer_thrs_%s_%s_70_minROI.mat',subject{ind},hemi{h},res{r}),'T1_D3_tmp');
%          save(sprintf('%s_nQSM_mface_mhand_mfoot_DDlayer_thrs_%s_%s_70_minROI.mat',subject{ind},hemi{h},res{r}),'T1_all_tmp');
% % 
%          % replace zero with NaN
% 
%          T1_D1_21(T1_D1_21==0) = NaN;
%          T1_D2_21(T1_D2_21==0) = NaN;
%          T1_D3_21(T1_D3_21==0) = NaN;
%          T1_all_21(T1_all_21==0) = NaN;
%          % 
%          T1_D1_21_neg(T1_D1_21_neg==0) = NaN;
%          T1_D2_21_neg(T1_D2_21_neg==0) = NaN;
%          T1_D3_21_neg(T1_D3_21_neg==0) = NaN;
%          T1_all_21_neg(T1_all_21_neg==0) = NaN;
% 
%          T1_D1_21_abs(T1_D1_21_abs==0) = NaN;
%          T1_D2_21_abs(T1_D2_21_abs==0) = NaN;
%          T1_D3_21_abs(T1_D3_21_abs==0) = NaN;
%          T1_all_21_abs(T1_all_21_abs==0) = NaN;
% 
%          T1_D1(T1_D1==0) = NaN;
%          T1_D2(T1_D2==0) = NaN;
%          T1_D3(T1_D3==0) = NaN;
%          T1_all(T1_all==0) = NaN;
% 
%          % calculate mean across vertices per layer (one value per layer)
%          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%          [rownum,colnum]=size(T1_D1_21);
% 
%          for ind5 = 1:rownum % iterate over all layers
%             mean_T1_D1_21(ind5,1) = nanmean(T1_D1_21(ind5,1:end)); %generates column vector with 4 rows (one row per layer containing average T1 value in the first column
%             mean_T1_D1_21(ind5,2) = 1; % generates columns vector with 4 rows containing finger number
%             mean_T1_D1_21(ind5,3) = ind5; % generates column vector with 4 rows containing layer number
%             mean_T1_D1_21(ind5,4) = subject_code(ind);
%             mean_T1_D1_21(ind5,5) = ind;
%             mean_T1_D1_21(ind5,6) = vertices_d1_S1;
%             mean_T1_D1_21(ind5,7) = h;
%             mean_T1_D1_21(ind5,8) = r;
% 
%             mean_T1_D2_21(ind5,1) = nanmean(T1_D2_21(ind5,1:end));
%             mean_T1_D2_21(ind5,2) = 2;
%             mean_T1_D2_21(ind5,3) = ind5;
%             mean_T1_D2_21(ind5,4) = subject_code(ind);
%             mean_T1_D2_21(ind5,5) = ind;
%             mean_T1_D2_21(ind5,6) = vertices_d2_S1;
%             mean_T1_D2_21(ind5,7) = h;
%             mean_T1_D2_21(ind5,8) = r;
% 
%             mean_T1_D3_21(ind5,1) = nanmean(T1_D3_21(ind5,1:end));
%             mean_T1_D3_21(ind5,2) = 3;
%             mean_T1_D3_21(ind5,3) = ind5;
%             mean_T1_D3_21(ind5,4) = subject_code(ind);
%             mean_T1_D3_21(ind5,5) = ind;
%             mean_T1_D3_21(ind5,6) = vertices_d3_S1;
%             mean_T1_D3_21(ind5,7) = h;
%             mean_T1_D3_21(ind5,8) = r;
% 
%             mean_T1_all_21(ind5,1) = nanmean(T1_all_21(ind5,1:end));
%             mean_T1_all_21(ind5,2) = 0;
%             mean_T1_all_21(ind5,3) = ind5;
%             mean_T1_all_21(ind5,4) = subject_code(ind);
%             mean_T1_all_21(ind5,5) = ind;
%             mean_T1_all_21(ind5,6) = vertices_d1_S1 + vertices_d2_S1 + vertices_d3_S1; %+vertices_d3_S1 + vertices_d4_S1 + vertices_d5_S1;
%             mean_T1_all_21(ind5,7) = h;
%             mean_T1_all_21(ind5,8) = r;
% 
%             mean_T1_D1_21_neg(ind5,1) = nanmean(T1_D1_21_neg(ind5,1:end)); %generates column vector with 4 rows (one row per layer containing average T1 value in the first column
%             mean_T1_D1_21_neg(ind5,2) = 1; % generates columns vector with 4 rows containing finger number
%             mean_T1_D1_21_neg(ind5,3) = ind5; % generates column vector with 4 rows containing layer number
%             mean_T1_D1_21_neg(ind5,4) = subject_code(ind);
%             mean_T1_D1_21_neg(ind5,5) = ind;
%             mean_T1_D1_21_neg(ind5,6) = vertices_d1_S1_neg;
%             mean_T1_D1_21_neg(ind5,7) = h;
%             mean_T1_D1_21_neg(ind5,8) = r;
% 
%             mean_T1_D2_21_neg(ind5,1) = nanmean(T1_D2_21_neg(ind5,1:end));
%             mean_T1_D2_21_neg(ind5,2) = 2;
%             mean_T1_D2_21_neg(ind5,3) = ind5;
%             mean_T1_D2_21_neg(ind5,4) = subject_code(ind);
%             mean_T1_D2_21_neg(ind5,5) = ind;
%             mean_T1_D2_21_neg(ind5,6) = vertices_d2_S1_neg;
%             mean_T1_D2_21_neg(ind5,7) = h;
%             mean_T1_D2_21_neg(ind5,8) = r;
% 
%             mean_T1_D3_21_neg(ind5,1) = nanmean(T1_D3_21_neg(ind5,1:end));
%             mean_T1_D3_21_neg(ind5,2) = 3;
%             mean_T1_D3_21_neg(ind5,3) = ind5;
%             mean_T1_D3_21_neg(ind5,4) = subject_code(ind);
%             mean_T1_D3_21_neg(ind5,5) = ind;
%             mean_T1_D3_21_neg(ind5,6) = vertices_d3_S1_neg;
%             mean_T1_D3_21_neg(ind5,7) = h;
%             mean_T1_D3_21_neg(ind5,8) = r;
% 
%             mean_T1_all_21_neg(ind5,1) = nanmean(T1_all_21_neg(ind5,1:end));
%             mean_T1_all_21_neg(ind5,2) = 0;
%             mean_T1_all_21_neg(ind5,3) = ind5;
%             mean_T1_all_21_neg(ind5,4) = subject_code(ind);
%             mean_T1_all_21_neg(ind5,5) = ind;
%             mean_T1_all_21_neg(ind5,6) = vertices_d1_S1_neg + vertices_d2_S1_neg + vertices_d3_S1_neg; %+vertices_d3_S1 + vertices_d4_S1 + vertices_d5_S1;
%             mean_T1_all_21_neg(ind5,7) = h;
%             mean_T1_all_21_neg(ind5,8) = r;
% 
%             mean_T1_D1_21_abs(ind5,1) = nanmean(T1_D1_21_abs(ind5,1:end)); %generates column vector with 4 rows (one row per layer containing average T1 value in the first column
%             mean_T1_D1_21_abs(ind5,2) = 1; % generates columns vector with 4 rows containing finger number
%             mean_T1_D1_21_abs(ind5,3) = ind5; % generates column vector with 4 rows containing layer number
%             mean_T1_D1_21_abs(ind5,4) = subject_code(ind);
%             mean_T1_D1_21_abs(ind5,5) = ind;
%             mean_T1_D1_21_abs(ind5,6) = vertices_d1_S1_abs;
%             mean_T1_D1_21_abs(ind5,7) = h;
%             mean_T1_D1_21_abs(ind5,8) = r;
% 
%             mean_T1_D2_21_abs(ind5,1) = nanmean(T1_D2_21_abs(ind5,1:end));
%             mean_T1_D2_21_abs(ind5,2) = 2;
%             mean_T1_D2_21_abs(ind5,3) = ind5;
%             mean_T1_D2_21_abs(ind5,4) = subject_code(ind);
%             mean_T1_D2_21_abs(ind5,5) = ind;
%             mean_T1_D2_21_abs(ind5,6) = vertices_d2_S1_abs;
%             mean_T1_D2_21_abs(ind5,7) = h;
%             mean_T1_D2_21_abs(ind5,8) = r;
% 
%             mean_T1_D3_21_abs(ind5,1) = nanmean(T1_D3_21_abs(ind5,1:end));
%             mean_T1_D3_21_abs(ind5,2) = 3;
%             mean_T1_D3_21_abs(ind5,3) = ind5;
%             mean_T1_D3_21_abs(ind5,4) = subject_code(ind);
%             mean_T1_D3_21_abs(ind5,5) = ind;
%             mean_T1_D3_21_abs(ind5,6) = vertices_d3_S1_abs;
%             mean_T1_D3_21_abs(ind5,7) = h;
%             mean_T1_D3_21_abs(ind5,8) = r;
% 
%             mean_T1_all_21_abs(ind5,1) = nanmean(T1_all_21_abs(ind5,1:end));
%             mean_T1_all_21_abs(ind5,2) = 0;
%             mean_T1_all_21_abs(ind5,3) = ind5;
%             mean_T1_all_21_abs(ind5,4) = subject_code(ind);
%             mean_T1_all_21_abs(ind5,5) = ind;
%             mean_T1_all_21_abs(ind5,6) = vertices_d1_S1_abs + vertices_d2_S1_abs + vertices_d3_S1_abs; %+vertices_d3_S1 + vertices_d4_S1 + vertices_d5_S1;
%             mean_T1_all_21_abs(ind5,7) = h;
%             mean_T1_all_21_abs(ind5,8) = r;
% 
%          end
% 
%          % write data to matfile
%          save(sprintf('%s_pQSM_mhand_21layer_mean_thrs_%s_%s_70_minROI.mat',subject{ind},hemi{h},res{r}),'mean_T1_D1_21');
%          save(sprintf('%s_pQSM_mface_21layer_mean_thrs_%s_%s_70_minROI.mat',subject{ind},hemi{h},res{r}),'mean_T1_D2_21');
%          save(sprintf('%s_pQSM_mfoot_21layer_mean_thrs_%s_%s_70_minROI.mat',subject{ind},hemi{h},res{r}),'mean_T1_D3_21');
%          save(sprintf('%s_pQSM_mface_mhand_mfoot_21layer_mean_thrs_%s_%s_70_minROI.mat',subject{ind},hemi{h},res{r}),'mean_T1_all_21');
% 
%          save(sprintf('%s_nQSM_mhand_21layer_mean_thrs_%s_%s_70_minROI.mat',subject{ind},hemi{h},res{r}),'mean_T1_D1_21_neg');
%          save(sprintf('%s_nQSM_mface_21layer_mean_thrs_%s_%s_70_minROI.mat',subject{ind},hemi{h},res{r}),'mean_T1_D2_21_neg');
%          save(sprintf('%s_nQSM_mfoot_21layer_mean_thrs_%s_%s_70_minROI.mat',subject{ind},hemi{h},res{r}),'mean_T1_D3_21_neg');
%          save(sprintf('%s_nQSM_mface_mhand_mfoot_21layer_mean_thrs_%s_%s_70_minROI.mat',subject{ind},hemi{h},res{r}),'mean_T1_all_21_neg');
% 
%          save(sprintf('%s_aQSM_mhand_21layer_mean_thrs_%s_%s_70_minROI.mat',subject{ind},hemi{h},res{r}),'mean_T1_D1_21_abs');
%          save(sprintf('%s_aQSM_mface_21layer_mean_thrs_%s_%s_70_minROI.mat',subject{ind},hemi{h},res{r}),'mean_T1_D2_21_abs');
%          save(sprintf('%s_aQSM_mfoot_21layer_mean_thrs_%s_%s_70_minROI.mat',subject{ind},hemi{h},res{r}),'mean_T1_D3_21_abs');
%          save(sprintf('%s_aQSM_mface_mhand_mfoot_21layer_mean_thrs_%s_%s_70_minROI.mat',subject{ind},hemi{h},res{r}),'mean_T1_all_21_abs');
% 
%          table_allsubjects_fingers_long_21 = [table_allsubjects_fingers_long_21; mean_T1_D1_21(:,:);mean_T1_D2_21(:,:);mean_T1_D3_21(:,:)];
%          table_allsubjects_fingers_long_21_neg = [table_allsubjects_fingers_long_21_neg; mean_T1_D1_21_neg(:,:);mean_T1_D2_21_neg(:,:);mean_T1_D3_21_neg(:,:)];
%          table_allsubjects_fingers_long_21_abs = [table_allsubjects_fingers_long_21_abs; mean_T1_D1_21_abs(:,:);mean_T1_D2_21_abs(:,:);mean_T1_D3_21_abs(:,:)];
% 
%          if h==1
%             table_allsubjects_fingers_short_21(ind,:) = [transpose(mean_T1_D1_21(:,1)) transpose(mean_T1_D2_21(:,1)) transpose(mean_T1_D3_21(:,1)) subject_code(ind) ind nanmean(mean_T1_D1_21(:,6)) nanmean(mean_T1_D2_21(:,6)) nanmean(mean_T1_D3_21(:,6)) h];
%             table_allsubjects_fingers_short_21_neg(ind,:) = [transpose(mean_T1_D1_21_neg(:,1)) transpose(mean_T1_D2_21_neg(:,1)) transpose(mean_T1_D3_21_neg(:,1)) subject_code(ind) ind nanmean(mean_T1_D1_21_neg(:,6)) nanmean(mean_T1_D2_21_neg(:,6)) nanmean(mean_T1_D3_21_neg(:,6)) h];
%             table_allsubjects_fingers_short_21_abs(ind,:) = [transpose(mean_T1_D1_21_abs(:,1)) transpose(mean_T1_D2_21_abs(:,1)) transpose(mean_T1_D3_21_abs(:,1)) subject_code(ind) ind nanmean(mean_T1_D1_21_abs(:,6)) nanmean(mean_T1_D2_21_abs(:,6)) nanmean(mean_T1_D3_21_abs(:,6)) h];
%          else
%             table_allsubjects_fingers_short_21_na(ind,:) = [transpose(mean_T1_D1_21(:,1)) transpose(mean_T1_D2_21(:,1)) transpose(mean_T1_D3_21(:,1)) subject_code(ind) ind nanmean(mean_T1_D1_21(:,6)) nanmean(mean_T1_D2_21(:,6)) nanmean(mean_T1_D3_21(:,6)) h];
%             table_allsubjects_fingers_short_21_neg_na(ind,:) = [transpose(mean_T1_D1_21_neg(:,1)) transpose(mean_T1_D2_21_neg(:,1)) transpose(mean_T1_D3_21_neg(:,1)) subject_code(ind) ind nanmean(mean_T1_D1_21_neg(:,6)) nanmean(mean_T1_D2_21_neg(:,6)) nanmean(mean_T1_D3_21_neg(:,6)) h];
%             table_allsubjects_fingers_short_21_abs_na(ind,:) = [transpose(mean_T1_D1_21_abs(:,1)) transpose(mean_T1_D2_21_abs(:,1)) transpose(mean_T1_D3_21_abs(:,1)) subject_code(ind) ind nanmean(mean_T1_D1_21_abs(:,6)) nanmean(mean_T1_D2_21_abs(:,6)) nanmean(mean_T1_D3_21_abs(:,6)) h];
% 
%          end
%          [rownum,colnum]=size(T1_D1);
%          for ind5 = 1:(rownum/3) % iterate over all layers
%             mean_T1_D1_abs(ind5,1) = nanmean(T1_D1(ind5,1:end)); %generates column vector with 4 rows (one row per layer containing average T1 value in the first column
%             mean_T1_D1_abs(ind5,2) = 1; % generates columns vector with 4 rows containing finger number
%             mean_T1_D1_abs(ind5,3) = ind5; % generates column vector with 4 rows containing layer number
%             mean_T1_D1_abs(ind5,4) = subject_code(ind);
%             mean_T1_D1_abs(ind5,5) = ind;
%             mean_T1_D1_abs(ind5,6) = vertices_d1_S1_abs;
%             mean_T1_D1_abs(ind5,7) = h;
% 
%             mean_T1_D2_abs(ind5,1) = nanmean(T1_D2(ind5,1:end));
%             mean_T1_D2_abs(ind5,2) = 2;
%             mean_T1_D2_abs(ind5,3) = ind5;
%             mean_T1_D2_abs(ind5,4) = subject_code(ind);
%             mean_T1_D2_abs(ind5,5) = ind;
%             mean_T1_D2_abs(ind5,6) = vertices_d2_S1_abs;
%             mean_T1_D2_abs(ind5,7) = h;
% 
%             mean_T1_D3_abs(ind5,1) = nanmean(T1_D3(ind5,1:end));
%             mean_T1_D3_abs(ind5,2) = 3;
%             mean_T1_D3_abs(ind5,3) = ind5;
%             mean_T1_D3_abs(ind5,4) = subject_code(ind);
%             mean_T1_D3_abs(ind5,5) = ind;
%             mean_T1_D3_abs(ind5,6) = vertices_d3_S1_abs;
%             mean_T1_D3_abs(ind5,7) = h;
% 
%             mean_T1_all_abs(ind5,1) = nanmean(T1_all(ind5,1:end));
%             mean_T1_all_abs(ind5,2) = 0;
%             mean_T1_all_abs(ind5,3) = ind5;
%             mean_T1_all_abs(ind5,4) = subject_code(ind);
%             mean_T1_all_abs(ind5,5) = ind;
%             mean_T1_all_abs(ind5,6) = vertices_d1_S1_abs + vertices_d2_S1_abs + vertices_d3_S1_abs; %+vertices_d3_S1 + vertices_d4_S1 + vertices_d5_S1;
%             mean_T1_all_abs(ind5,7) = h;
% 
%             mean_T1_D1_pos(ind5,1) = nanmean(T1_D1(ind5+3,1:end)); %generates column vector with 4 rows (one row per layer containing average T1 value in the first column
%             mean_T1_D1_pos(ind5,2) = 1; % generates columns vector with 4 rows containing finger number
%             mean_T1_D1_pos(ind5,3) = ind5; % generates column vector with 4 rows containing layer number
%             mean_T1_D1_pos(ind5,4) = subject_code(ind);
%             mean_T1_D1_pos(ind5,5) = ind;
%             mean_T1_D1_pos(ind5,6) = vertices_d1_S1;
%             mean_T1_D1_pos(ind5,7) = h;
% 
%             mean_T1_D2_pos(ind5,1) = nanmean(T1_D2(ind5+3,1:end));
%             mean_T1_D2_pos(ind5,2) = 2;
%             mean_T1_D2_pos(ind5,3) = ind5;
%             mean_T1_D2_pos(ind5,4) = subject_code(ind);
%             mean_T1_D2_pos(ind5,5) = ind;
%             mean_T1_D2_pos(ind5,6) = vertices_d2_S1;
%             mean_T1_D2_pos(ind5,7) = h;
% 
%             mean_T1_D3_pos(ind5,1) = nanmean(T1_D3(ind5+3,1:end));
%             mean_T1_D3_pos(ind5,2) = 3;
%             mean_T1_D3_pos(ind5,3) = ind5;
%             mean_T1_D3_pos(ind5,4) = subject_code(ind);
%             mean_T1_D3_pos(ind5,5) = ind;
%             mean_T1_D3_pos(ind5,6) = vertices_d3_S1;
%             mean_T1_D3_pos(ind5,7) = h;
% 
%             mean_T1_all_pos(ind5,1) = nanmean(T1_all(ind5+3,1:end));
%             mean_T1_all_pos(ind5,2) = 0;
%             mean_T1_all_pos(ind5,3) = ind5;
%             mean_T1_all_pos(ind5,4) = subject_code(ind);
%             mean_T1_all_pos(ind5,5) = ind;
%             mean_T1_all_pos(ind5,6) = vertices_d1_S1 + vertices_d2_S1 + vertices_d3_S1; %+vertices_d3_S1 + vertices_d4_S1 + vertices_d5_S1;
%             mean_T1_all_pos(ind5,7) = h;
% 
%             mean_T1_D1_neg(ind5,1) = nanmean(T1_D1(ind5+6,1:end)); %generates column vector with 4 rows (one row per layer containing average T1 value in the first column
%             mean_T1_D1_neg(ind5,2) = 1; % generates columns vector with 4 rows containing finger number
%             mean_T1_D1_neg(ind5,3) = ind5; % generates column vector with 4 rows containing layer number
%             mean_T1_D1_neg(ind5,4) = subject_code(ind);
%             mean_T1_D1_neg(ind5,5) = ind;
%             mean_T1_D1_neg(ind5,6) = vertices_d1_S1_neg;
%             mean_T1_D1_neg(ind5,7) = h;
% 
%             mean_T1_D2_neg(ind5,1) = nanmean(T1_D2(ind5+6,1:end));
%             mean_T1_D2_neg(ind5,2) = 2;
%             mean_T1_D2_neg(ind5,3) = ind5;
%             mean_T1_D2_neg(ind5,4) = subject_code(ind);
%             mean_T1_D2_neg(ind5,5) = ind;
%             mean_T1_D2_neg(ind5,6) = vertices_d2_S1_neg;
%             mean_T1_D2_neg(ind5,7) = h;
% 
%             mean_T1_D3_neg(ind5,1) = nanmean(T1_D3(ind5+6,1:end));
%             mean_T1_D3_neg(ind5,2) = 3;
%             mean_T1_D3_neg(ind5,3) = ind5;
%             mean_T1_D3_neg(ind5,4) = subject_code(ind);
%             mean_T1_D3_neg(ind5,5) = ind;
%             mean_T1_D3_neg(ind5,6) = vertices_d3_S1_neg;
%             mean_T1_D3_neg(ind5,7) = h;
% 
%             mean_T1_all_neg(ind5,1) = nanmean(T1_all(ind5+6,1:end));
%             mean_T1_all_neg(ind5,2) = 0;
%             mean_T1_all_neg(ind5,3) = ind5;
%             mean_T1_all_neg(ind5,4) = subject_code(ind);
%             mean_T1_all_neg(ind5,5) = ind;
%             mean_T1_all_neg(ind5,6) = vertices_d1_S1_neg + vertices_d2_S1_neg + vertices_d3_S1_neg; %+vertices_d3_S1 + vertices_d4_S1 + vertices_d5_S1;
%             mean_T1_all_neg(ind5,7) = h;
%          end
% 
%          table_allsubjects_fingers_long_abs = [table_allsubjects_fingers_long_abs; mean_T1_D1_abs(:,:);mean_T1_D2_abs(:,:);mean_T1_D3_abs(:,:)];
%          table_allsubjects_fingers_long_pos = [table_allsubjects_fingers_long_pos; mean_T1_D1_pos(:,:);mean_T1_D2_pos(:,:);mean_T1_D3_pos(:,:)];
%          table_allsubjects_fingers_long_neg = [table_allsubjects_fingers_long_neg; mean_T1_D1_neg(:,:);mean_T1_D2_neg(:,:);mean_T1_D3_neg(:,:)];
% 
%          %l=length(table_allsubjects_fingers_short(ind,:))
%          if h==1
%             table_allsubjects_fingers_short_abs(ind,:) = [transpose(mean_T1_D1_abs(:,1)) transpose(mean_T1_D2_abs(:,1)) transpose(mean_T1_D3_abs(:,1)) subject_code(ind) ind nanmean(mean_T1_D1_abs(:,6)) nanmean(mean_T1_D2_abs(:,6)) nanmean(mean_T1_D3_abs(:,6)) h];
%             table_allsubjects_fingers_short_pos(ind,:) = [transpose(mean_T1_D1_pos(:,1)) transpose(mean_T1_D2_pos(:,1)) transpose(mean_T1_D3_pos(:,1)) subject_code(ind) ind nanmean(mean_T1_D1_pos(:,6)) nanmean(mean_T1_D2_pos(:,6)) nanmean(mean_T1_D3_pos(:,6)) h];
%             table_allsubjects_fingers_short_neg(ind,:) = [transpose(mean_T1_D1_neg(:,1)) transpose(mean_T1_D2_neg(:,1)) transpose(mean_T1_D3_neg(:,1)) subject_code(ind) ind nanmean(mean_T1_D1_neg(:,6)) nanmean(mean_T1_D2_neg(:,6)) nanmean(mean_T1_D3_neg(:,6)) h];
% 
%          else
%             table_allsubjects_fingers_short_na_abs(ind,:) = [transpose(mean_T1_D1_abs(:,1)) transpose(mean_T1_D2_abs(:,1)) transpose(mean_T1_D3_abs(:,1)) subject_code(ind) ind nanmean(mean_T1_D1_abs(:,6)) nanmean(mean_T1_D2_abs(:,6)) nanmean(mean_T1_D3_abs(:,6)) h];
%             table_allsubjects_fingers_short_na_pos(ind,:) = [transpose(mean_T1_D1_pos(:,1)) transpose(mean_T1_D2_pos(:,1)) transpose(mean_T1_D3_pos(:,1)) subject_code(ind) ind nanmean(mean_T1_D1_pos(:,6)) nanmean(mean_T1_D2_pos(:,6)) nanmean(mean_T1_D3_pos(:,6)) h];
%             table_allsubjects_fingers_short_na_neg(ind,:) = [transpose(mean_T1_D1_neg(:,1)) transpose(mean_T1_D2_neg(:,1)) transpose(mean_T1_D3_neg(:,1)) subject_code(ind) ind nanmean(mean_T1_D1_neg(:,6)) nanmean(mean_T1_D2_neg(:,6)) nanmean(mean_T1_D3_neg(:,6)) h];
% 
%          end
% % % 
%          save(sprintf('%s_aQSM_mhand_DDlayer_mean_thrs_%s_%s_70_minROI.mat',subject{ind},hemi{h},res{r}),'mean_T1_D1_abs');
%          save(sprintf('%s_aQSM_mface_DDlayer_mean_thrs_%s_%s_70_minROI.mat',subject{ind},hemi{h},res{r}),'mean_T1_D2_abs');
%          save(sprintf('%s_aQSM_mfoot_DDlayer_mean_thrs_%s_%s_70_minROI.mat',subject{ind},hemi{h},res{r}),'mean_T1_D3_abs');
%          save(sprintf('%s_aQSM_mface_mhand_mfoot_DDlayer_mean_thrs_%s_%s_70_minROI.mat',subject{ind},hemi{h},res{r}),'mean_T1_all_abs');
% 
%          save(sprintf('%s_pQSM_mhand_DDlayer_mean_thrs_%s_%s_70_minROI.mat',subject{ind},hemi{h},res{r}),'mean_T1_D1_pos');
%          save(sprintf('%s_pQSM_mface_DDlayer_mean_thrs_%s_%s_70_minROI.mat',subject{ind},hemi{h},res{r}),'mean_T1_D2_pos');
%          save(sprintf('%s_pQSM_mfoot_DDlayer_mean_thrs_%s_%s_70_minROI.mat',subject{ind},hemi{h},res{r}),'mean_T1_D3_pos');
%          save(sprintf('%s_pQSM_mface_mhand_mfoot_DDlayer_mean_thrs_%s_%s_70_minROI.mat',subject{ind},hemi{h},res{r}),'mean_T1_all_pos');
% 
%          save(sprintf('%s_nQSM_mhand_DDlayer_mean_thrs_%s_%s_70_minROI.mat',subject{ind},hemi{h},res{r}),'mean_T1_D1_neg');
%          save(sprintf('%s_nQSM_mface_DDlayer_mean_thrs_%s_%s_70_minROI.mat',subject{ind},hemi{h},res{r}),'mean_T1_D2_neg');
%          save(sprintf('%s_nQSM_mfoot_DDlayer_mean_thrs_%s_%s_70_minROI.mat',subject{ind},hemi{h},res{r}),'mean_T1_D3_neg');
%          save(sprintf('%s_nQSM_mface_mhand_mfoot_DDlayer_mean_thrs_%s_%s_70_minROI.mat',subject{ind},hemi{h},res{r}),'mean_T1_all_neg');
        
     % else
     %    warning('Nothing done. Getting next subject...');
     % end   

    end

    end

end
% C_abs = [table_allsubjects_fingers_short_abs table_allsubjects_fingers_short_na_abs]
% C_pos = [table_allsubjects_fingers_short_pos table_allsubjects_fingers_short_na_pos]
% C_neg = [table_allsubjects_fingers_short_neg table_allsubjects_fingers_short_na_neg]
% 
% C_abs_21 = [table_allsubjects_fingers_short_21_abs table_allsubjects_fingers_short_21_abs_na]
% C_pos_21 = [table_allsubjects_fingers_short_21 table_allsubjects_fingers_short_21_na]
% C_neg_21 = [table_allsubjects_fingers_short_21_neg table_allsubjects_fingers_short_21_neg_na]
% 
% save(sprintf('allsubjects_pQSM_mFaceHandFoot_short_21Layers_%s_70_minROI.mat',res{1}),'C_pos_21');
% save(sprintf('allsubjects_pQSM_mFaceHandFoot_long_21Layers_%s_70_minROI.mat',res{1}),'table_allsubjects_fingers_long_21');
% csvwrite(sprintf('allsubjects_pQSM_mFaceHandFoot_short_21Layers_%s_70_minROI.csv',res{1}),C_pos_21);
% csvwrite(sprintf('allsubjects_pQSM_mFaceHandFoot_long_21Layers_%s_70_minROI.csv',res{1}),table_allsubjects_fingers_long_21);
% 
% save(sprintf('allsubjects_nQSM_mFaceHandFoot_short_21Layers_%s_70_minROI.mat',res{1}),'C_neg_21');
% save(sprintf('allsubjects_nQSM_mFaceHandFoot_long_21Layers_%s_70_minROI.mat',res{1}),'table_allsubjects_fingers_long_21_neg');
% csvwrite(sprintf('allsubjects_nQSM_mFaceHandFoot_short_21Layers_%s_70_minROI.csv',res{1}),C_neg_21);
% csvwrite(sprintf('allsubjects_nQSM_mFaceHandFoot_long_21Layers_%s_70_minROI.csv',res{1}),table_allsubjects_fingers_long_21_neg);
% 
% save(sprintf('allsubjects_aQSM_mFaceHandFoot_short_21Layers_%s_70_minROI.mat',res{1}),'C_abs_21');
% save(sprintf('allsubjects_aQSM_mFaceHandFoot_long_21Layers_%s_70_minROI.mat',res{1}),'table_allsubjects_fingers_long_21_abs');
% csvwrite(sprintf('allsubjects_aQSM_mFaceHandFoot_short_21Layers_%s_70_minROI.csv',res{1}),C_abs_21);
% csvwrite(sprintf('allsubjects_aQSM_mFaceHandFoot_long_21Layers_%s_70_minROI.csv',res{1}),table_allsubjects_fingers_long_21_abs);
% 
% 
% 
% save(sprintf('allsubjects_pQSM_mFaceHandFoot_short_DDLayers_%s_70_minROI.mat',res{1}),'C_pos');
% save(sprintf('allsubjects_pQSM_mFaceHandFoot_long_DDLayers_%s_70_minROI.mat',res{1}),'table_allsubjects_fingers_long_pos');
% csvwrite(sprintf('allsubjects_pQSM_mFaceHandFoot_short_DDLayers_%s_70_minROI.csv',res{1}),C_pos);
% csvwrite(sprintf('allsubjects_pQSM_mFaceHandFoot_long_DDLayers_%s_70_minROI.csv',res{1}),table_allsubjects_fingers_long_pos);
% 
% save(sprintf('allsubjects_nQSM_mFaceHandFoot_short_DDLayers_%s_70_minROI.mat',res{1}),'C_neg');
% save(sprintf('allsubjects_nQSM_mFaceHandFoot_long_DDLayers_%s_70_minROI.mat',res{1}),'table_allsubjects_fingers_long_neg');
% csvwrite(sprintf('allsubjects_nQSM_mFaceHandFoot_short_DDLayers_%s_70_minROI.csv',res{1}),C_neg);
% csvwrite(sprintf('allsubjects_nQSM_mFaceHandFoot_long_DDLayers_%s_70_minROI.csv',res{1}),table_allsubjects_fingers_long_neg);
% 
% save(sprintf('allsubjects_aQSM_mFaceHandFoot_short_DDLayers_%s_70_minROI.mat',res{1}),'C_abs');
% save(sprintf('allsubjects_aQSM_mFaceHandFoot_long_DDLayers_%s_70_minROI.mat',res{1}),'table_allsubjects_fingers_long_abs');
% csvwrite(sprintf('allsubjects_aQSM_mFaceHandFoot_short_DDLayers_%s_70_minROI.csv',res{1}),C_abs);
% csvwrite(sprintf('allsubjects_aQSM_mFaceHandFoot_long_DDLayers_%s_70_minROI.csv',res{1}),table_allsubjects_fingers_long_abs);
