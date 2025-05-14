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
% T1file_nameSurface = '*_surf_2-2_inf__3-6.vtk';
headerdir = '/home/juli/projects/onehander/data';

function_nameSurface = {'*_hand_S1_smoothdata.vtk','*_face_S1_smoothdata.vtk','*_foot_S1_smoothdata.vtk'};

% definitions
% for laptop
datadir = '/home/juli/projects/onehander/data';
datadir_func = '/home/juli/projects/onehander/data';

script_dir = '/home/juli/projects/onehander/scripts';
result_dir = '/home/juli/projects/onehander/results/ReadOutLayerData_T1_HandFaceFoot';

layer_def = '/home/juli/projects/onehander/results/extract_Anatomical_layers';

addpath(script_dir)

% specifiy number of functional maps
n_fingers = 3;

file_ext = {'qT1_S1_inner_layers','qT1_S1_middle_layers','qT1_S1_outer_layers'};
%file_ext = {'qT1_S1_inner_ylayers','qT1_S1_middle_ylayers','qT1_S1_outer_ylayers'};

file_ext_loc = {'hand','face','foot'};

hemi = {'left','right'};

%res = {'orig','resampled'};
res = {'orig'};

table_allsubjects_fingers_long = [];
table_allsubjects_fingers_short = [];

table_allsubjects_fingers_long_21 = [];
table_allsubjects_fingers_short_21 = [];

%% start loop to: extract T1 values and t-values for all subjects (one iteration per subject subject), store and plot data
for ind = 2:length(subject)

    for h = 1:length(hemi)  

        if ind == 2 && h == 1
            function_nameSurfaceFolder = {'*ECAA','*EDAA','*EHAA'};
        else
            function_nameSurfaceFolder = {'*ECAA','*EDAA','*EEAA'};
        end        

        for r = 1:length(res)

            % if h ==1
            %     function_nameSurfaceFolder = sprintf('%s/%s/processed/output/Mapping_WB_S1_FINAL/exp-0000/exp-0000-AAAAAAAA/SurfaceMeshMapping',datadir,subject{ind});
            % else
            %     function_nameSurfaceFolder = sprintf('%s/%s/processed/output/Mapping_WB_S1_FINAL/exp-0000/exp-0000-ABAAAAAA/SurfaceMeshMapping',datadir,subject{ind});
            % end
        
            % cd(function_nameSurfaceFolder)
            % T1fileinfo =  dir(sprintf('**/%s',T1file_nameSurface));
            % T1file = fullfile(T1fileinfo(1).folder, T1fileinfo(1).name)
            % 
            % [vertex,face,qT1vals,header1,header2,header3] = read_vtk(T1file);
    
             clear('T1_layers_area3b')
             clear('T1_layers_area3b_bin')
             clear('func_tval_finger')
             clear('finger_winner')
             clear('finger_winner_S1')
        
             % load layer definition
             if h == 1
                 load(sprintf('%s/%s/idx_innerLayer_l_a_70.mat',layer_def,subject{ind}));
                 load(sprintf('%s/%s/idx_middleLayer_l_a_70.mat',layer_def,subject{ind}));
                 load(sprintf('%s/%s/idx_outerLayer_l_a_70.mat',layer_def,subject{ind}));
                 ind_inner = ind_inner;
                 ind_middle = ind_middle;
                 ind_outer = ind_outer;
             elseif h == 2
                 load(sprintf('%s/%s/idx_innerLayer_r_na_70.mat',layer_def,subject{ind}));
                 load(sprintf('%s/%s/idx_middleLayer_r_na_70.mat',layer_def,subject{ind}));
                 load(sprintf('%s/%s/idx_outerLayer_r_na_70.mat',layer_def,subject{ind}));
                 ind_inner = ind_inner; %ind_inner_o
                 ind_middle = ind_middle; %ind_middle_o
                 ind_outer = ind_outer; %ind_outer_o
             end
            
             funcfile_folder = 'exp-0000'
        
             % go to T1file folder and load extracted qT1 values for 21 layers
             cd (sprintf('%s/%s/processed/output/ReadOutLayerData_T1_21Layers',datadir, subject{ind}));
             load(sprintf('%s_T1_BA3b_21layer_%s_-2.8_%s.mat',subject{ind},hemi{h},res{r}));

             % area3b = T1_all;
             % 
             % T1_layers_area3b_bin = area3b(10,:);
             % T1_layers_area3b_bin(T1_layers_area3b_bin>0) = 1;
             % 
             % load(sprintf('%s_T1_WB_21layer_%s_-2.8.mat',subject{ind},hemi{h}));

             area3b = T1_all;
             T1_layers_area3b = [];

             % replace zero with NaN
             T1_layers_area3b_NaN = area3b;
             T1_layers_area3b_NaN(T1_layers_area3b_NaN==0) = NaN;
                   
             % average qT1 values into anatomically relevant layer compartments
             ind_inner(1,1)
             ind_inner(1,2)
             T1_layers_area3b(1,1:length(area3b)) = nanmean(T1_layers_area3b_NaN(ind_inner(1,1):ind_inner(1,2),:),1); %returns a row vector containing the mean for each column.
             ind_middle(1,1)
             ind_middle(1,2)
             T1_layers_area3b(2,1:length(area3b)) = nanmean(T1_layers_area3b_NaN(ind_middle(1,1):ind_middle(1,2),:),1); %returns a row vector containing the mean for each column.
             ind_outer(1,1)
             ind_outer(1,2)
             T1_layers_area3b(3,1:length(area3b)) = nanmean(T1_layers_area3b_NaN(ind_outer(1,1):ind_outer(1,2),:),1); %returns a row vector containing the mean for each column.
             
             % T1_layers_area3b(isnan(T1_layers_area3b)) = 0;
             % 
             % T1_layers_area3b_fin = T1_layers_area3b .* T1_layers_area3b_bin;
             % 
             % filename_new = sprintf('%s_qT1_S1_inner_DDlayers_%s_%s.vtk',subject{ind},hemi{h},res{r}); 
             % write_vtk(filename_new,header1,header2,header3,vertex,face,T1_layers_area3b_fin(1,:));
             % 
             % filename_new = sprintf('%s_qT1_S1_middle_DDlayers_%s_%s.vtk',subject{ind},hemi{h},res{r}); 
             % write_vtk(filename_new,header1,header2,header3,vertex,face,T1_layers_area3b_fin(2,:));
             % 
             % filename_new = sprintf('%s_qT1_S1_outer_DDlayers_%s_%s.vtk',subject{ind},hemi{h},res{r}); 
             % write_vtk(filename_new,header1,header2,header3,vertex,face,T1_layers_area3b_fin(3,:));


     % Calculate mean value per layer compartment     
     mean_T1_layers_3b(ind,1) = nanmean(T1_layers_area3b(1,:));
     mean_T1_layers_3b(ind,2) = nanmean(T1_layers_area3b(2,:));
     mean_T1_layers_3b(ind,3) = nanmean(T1_layers_area3b(3,:));
     % 
     % % generate binary S1 mask
     T1_layers_area3b_bin = area3b(10,:);
     T1_layers_area3b_bin(isnan(T1_layers_area3b_bin)) = 0;
     T1_layers_area3b_bin(T1_layers_area3b_bin>0) = 1;

     % start loop to extract values from functional maps: read out prfCL values for all fingers (one iteration per finger)
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for ind3 = 1:n_fingers

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

     if face_hand == 1

         % for win = 1:size(func_tval_finger,1)
         %     tmp_sort = sort(func_tval_finger(win,:));
         %     top_3000 = tmp_sort(end-2999)
         % 
         %     tmp_win = func_tval_finger(win,:);
         %     tmp_win(tmp_win<top_3000) = 0;
         % 
         %     tmp = tmp_win;
         %     tmp = tmp(tmp>0);
         %     length(tmp)
         % 
         %     func_tval_finger(win,:) = tmp_win;
         % end

         % winner takes it all: takes care of overlapping areas (introduced by remapping divided prf digit maps and smoothing), therefore only keep the maximum finger value if there are multiple values per voxel (multiple values per column), set all other finger values in the same row to zero 
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         % M = vector including maximum values of each column; I = vector including row indeces of maximum values
         % M also gives mask for all five fingers/hand
         [M,I] = max(func_tval_finger);
         M_bin = M;
         M_bin(M_bin~=0) = 1;

         % to separate fingers into different rows define winner array of zeros
         finger_winner = zeros(n_fingers,length(M));

         % start loop to fill matrix of zeros with maximum prf var values extacted from the different finger maps
         for ind5 = 1:length(M)
            finger_winner(I(ind5),ind5) = M(ind5);
         end

         length_winner = [];

         for win = 1:size(finger_winner,1)
             tmp = finger_winner(win,:);
             tmp = tmp(tmp>0);
             length_winner(win) = length(tmp);
         end

         for win = 1:size(finger_winner,1)
             tmp_sort = sort(finger_winner(win,:));
             %min(length_winner)
             top90 = round(length_winner(win)*0.7)
             %top = tmp_sort(end-1788)%tmp_sort(end-(min(top90)-1)) %min vertex size of 70% highest vertices identified before, script needs to bechanged to process hemis together if smallest ROI size is to be identified here
             top = tmp_sort(end-top90)

             tmp_win = finger_winner(win,:);
             tmp_win(tmp_win<top) = 0;

             tmp = tmp_win;
             tmp = tmp(tmp>0);
             length(tmp)

             finger_winner(win,:) = tmp_win;

             finger_winner_zvals = vertex(3,finger_winner(win,:)>0);

             %min(finger_winner_zvals)

             finger_winner_zvals_mean = mean(finger_winner_zvals);
             finger_winner_zvals_sd = std(finger_winner_zvals);

             %test = 'finger_winner_tmp'

             finger_winner_tmp = finger_winner(win,:);
             finger_winner_tmp(1,vertex(3,:)>finger_winner_zvals_mean+(2*finger_winner_zvals_sd)) = 0;
             finger_winner_tmp(1,vertex(3,:)<finger_winner_zvals_mean-(2*finger_winner_zvals_sd)) = 0;

             finger_winner(win,:) = finger_winner_tmp;

             %if ind3 == 3
                finger_winner_yvals = vertex(2,finger_winner(win,:)>0);
                finger_winner_yvals_mean = mean(finger_winner_yvals);
                finger_winner_yvals_sd = std(finger_winner_yvals);

                finger_winner_tmp(1,vertex(2,:)>finger_winner_yvals_mean+(2*finger_winner_yvals_sd)) = 0;
                finger_winner_tmp(1,vertex(2,:)<finger_winner_yvals_mean-(2*finger_winner_yvals_sd)) = 0;

                finger_winner(win,:) = finger_winner_tmp;

                finger_winner_xvals = vertex(1,finger_winner(win,:)>0);
                finger_winner_xvals_mean = mean(finger_winner_xvals);
                finger_winner_xvals_sd = std(finger_winner_xvals);

                finger_winner_tmp(1,vertex(1,:)>finger_winner_xvals_mean+(2*finger_winner_xvals_sd)) = 0;
                finger_winner_tmp(1,vertex(1,:)<finger_winner_xvals_mean-(2*finger_winner_xvals_sd)) = 0;

                finger_winner(win,:) = finger_winner_tmp;
             %end

             %min(vertex(2,finger_winner(win,:)>0))
         end

         cd (result_dir);

         %binarize
         finger_winner_bin = finger_winner;
         finger_winner_bin(finger_winner_bin~=0) = 1;

         sum(finger_winner_bin(1,:)>0)
         sum(finger_winner_bin(2,:)>0)
         sum(finger_winner_bin(3,:)>0)

         % apply S1 mask to final finger maps
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         finger_winner_S1 = zeros(n_fingers,length(finger_winner));

         % for i = 1:n_fingers
         %    finger_winner_S1(i,:) = finger_winner(i,:);
         % 
         %    % filename_new = sprintf('%s_%s_win_thrs_%s_%s.vtk',subject{ind},file_ext_loc{i},hemi{h},res{r}); 
         %    % write_vtk(filename_new,header1,header2,header3,vertex,face,finger_winner_S1(i,:))
         % end

         % write data to matfile
         save(sprintf('%s_body_part_win_S1_thrs_%s_%s_70.mat',subject{ind},hemi{h},res{r}),'finger_winner_S1');

         %binarize
         finger_winner_S1_bin = finger_winner_S1;
         finger_winner_S1_bin(finger_winner_S1_bin~=0) = 1;

         % multiply binarized finger masks with T1 values
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         T1_D1_21 = area3b.*finger_winner_bin(1,:); %d1 21 layers
         T1_D2_21 = area3b.*finger_winner_bin(2,:); %d2 21 layers
         T1_D3_21 = area3b.*finger_winner_bin(3,:); %d2 21 layers
         T1_all_21 = area3b.*M_bin; %hand (all fingers) %all 21 layers  

         T1_D1 = T1_layers_area3b.*finger_winner_bin(1,:); %d1 DD layers
         T1_D2 = T1_layers_area3b.*finger_winner_bin(2,:); %d2 DD layers
         T1_D3 = T1_layers_area3b.*finger_winner_bin(3,:); %d2 DD layers
         T1_all = T1_layers_area3b.*M_bin; %hand (all fingers) %all DD layers

         % calculate # of survived vertices for each finger map
         vertices_d1_S1 = sum(T1_D1_21(10,:)>0)
         vertices_d2_S1 = sum(T1_D2_21(10,:)>0)
         vertices_d3_S1 = sum(T1_D3_21(10,:)>0)

         % write data to matfile

         save(sprintf('%s_T1_mhand_21layer_thrs_%s_%s_70.mat',subject{ind},hemi{h},res{r}),'T1_D1_21');
         save(sprintf('%s_T1_mface_21layer_thrs_%s_%s_70.mat',subject{ind},hemi{h},res{r}),'T1_D2_21');
         save(sprintf('%s_T1_mfoot_21layer_thrs_%s_%s_70.mat',subject{ind},hemi{h},res{r}),'T1_D3_21');
         save(sprintf('%s_T1_mface_mhand_mfoot_21layer_thrs_%s_%s_70.mat',subject{ind},hemi{h},res{r}),'T1_all_21');
% 
         % save(sprintf('%s_T1_mhand_DDlayer_yLayers_thrs_70.mat',subject{ind}),'T1_D1');
         % save(sprintf('%s_T1_mface_DDlayer_yLayers_thrs_70.mat',subject{ind}),'T1_D2');
         % save(sprintf('%s_T1_mfoot_DDlayer_yLayers_thrs_70.mat',subject{ind}),'T1_D3');
         % save(sprintf('%s_T1_mface_mhand_mfoot_DDlayer_yLayers_thrs_70.mat',subject{ind}),'T1_all');

         save(sprintf('%s_T1_mhand_DDlayer_thrs_%s_%s_70.mat',subject{ind},hemi{h},res{r}),'T1_D1');
         save(sprintf('%s_T1_mface_DDlayer_thrs_%s_%s_70.mat',subject{ind},hemi{h},res{r}),'T1_D2');
         save(sprintf('%s_T1_mfoot_DDlayer_thrs_%s_%s_70.mat',subject{ind},hemi{h},res{r}),'T1_D3');
         save(sprintf('%s_T1_mface_mhand_mfoot_DDlayer_thrs_%s_%s_70.mat',subject{ind},hemi{h},res{r}),'T1_all');

         % replace zero with NaN

         T1_D1_21(T1_D1_21==0) = NaN;
         T1_D2_21(T1_D2_21==0) = NaN;
         T1_D3_21(T1_D3_21==0) = NaN;
         T1_all_21(T1_all_21==0) = NaN;
         % 
         T1_D1(T1_D1==0) = NaN;
         T1_D2(T1_D2==0) = NaN;
         T1_D3(T1_D3==0) = NaN;
         T1_all(T1_all==0) = NaN;

         % calculate mean across vertices per layer (one value per layer)
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         [rownum,colnum]=size(T1_D1_21);

         for ind5 = 1:rownum % iterate over all layers
            mean_T1_D1_21(ind5,1) = nanmean(T1_D1_21(ind5,1:end)); %generates column vector with 4 rows (one row per layer containing average T1 value in the first column
            mean_T1_D1_21(ind5,2) = 1; % generates columns vector with 4 rows containing finger number
            mean_T1_D1_21(ind5,3) = ind5; % generates column vector with 4 rows containing layer number
            mean_T1_D1_21(ind5,4) = subject_code(ind);
            mean_T1_D1_21(ind5,5) = ind;
            mean_T1_D1_21(ind5,6) = vertices_d1_S1;
            mean_T1_D1_21(ind5,7) = h;
            mean_T1_D1_21(ind5,8) = r;

            mean_T1_D2_21(ind5,1) = nanmean(T1_D2_21(ind5,1:end));
            mean_T1_D2_21(ind5,2) = 2;
            mean_T1_D2_21(ind5,3) = ind5;
            mean_T1_D2_21(ind5,4) = subject_code(ind);
            mean_T1_D2_21(ind5,5) = ind;
            mean_T1_D2_21(ind5,6) = vertices_d2_S1;
            mean_T1_D2_21(ind5,7) = h;
            mean_T1_D2_21(ind5,8) = r;

            mean_T1_D3_21(ind5,1) = nanmean(T1_D3_21(ind5,1:end));
            mean_T1_D3_21(ind5,2) = 3;
            mean_T1_D3_21(ind5,3) = ind5;
            mean_T1_D3_21(ind5,4) = subject_code(ind);
            mean_T1_D3_21(ind5,5) = ind;
            mean_T1_D3_21(ind5,6) = vertices_d3_S1;
            mean_T1_D3_21(ind5,7) = h;
            mean_T1_D3_21(ind5,8) = r;

            mean_T1_all_21(ind5,1) = nanmean(T1_all_21(ind5,1:end));
            mean_T1_all_21(ind5,2) = 0;
            mean_T1_all_21(ind5,3) = ind5;
            mean_T1_all_21(ind5,4) = subject_code(ind);
            mean_T1_all_21(ind5,5) = ind;
            mean_T1_all_21(ind5,6) = vertices_d1_S1 + vertices_d2_S1 + vertices_d3_S1; %+vertices_d3_S1 + vertices_d4_S1 + vertices_d5_S1;
            mean_T1_all_21(ind5,7) = h;
            mean_T1_all_21(ind5,8) = r;

         end

         % write data to matfile
         save(sprintf('%s_T1_mhand_21layer_mean_thrs_%s_%s_70.mat',subject{ind},hemi{h},res{r}),'mean_T1_D1_21');
         save(sprintf('%s_T1_mface_21layer_mean_thrs_%s_%s_70.mat',subject{ind},hemi{h},res{r}),'mean_T1_D2_21');
         save(sprintf('%s_T1_mfoot_21layer_mean_thrs_%s_%s_70.mat',subject{ind},hemi{h},res{r}),'mean_T1_D3_21');
         save(sprintf('%s_T1_mface_mhand_mfoot_21layer_mean_thrs_%s_%s_70.mat',subject{ind},hemi{h},res{r}),'mean_T1_all_21');

         table_allsubjects_fingers_long_21 = [table_allsubjects_fingers_long_21; mean_T1_D1_21(:,:);mean_T1_D2_21(:,:);mean_T1_D3_21(:,:)];

         if h==1
            table_allsubjects_fingers_short_21(ind,:) = [transpose(mean_T1_D1_21(:,1)) transpose(mean_T1_D2_21(:,1)) transpose(mean_T1_D3_21(:,1)) subject_code(ind) ind nanmean(mean_T1_D1_21(:,6)) nanmean(mean_T1_D2_21(:,6)) nanmean(mean_T1_D3_21(:,6)) h];
         else
            table_allsubjects_fingers_short_21_na(ind,:) = [transpose(mean_T1_D1_21(:,1)) transpose(mean_T1_D2_21(:,1)) transpose(mean_T1_D3_21(:,1)) subject_code(ind) ind nanmean(mean_T1_D1_21(:,6)) nanmean(mean_T1_D2_21(:,6)) nanmean(mean_T1_D3_21(:,6)) h];
         end

         [rownum,colnum]=size(T1_D1);

         for ind5 = 1:rownum % iterate over all layers
            mean_T1_D1(ind5,1) = nanmean(T1_D1(ind5,1:end)); %generates column vector with 4 rows (one row per layer containing average T1 value in the first column
            mean_T1_D1(ind5,2) = 1; % generates columns vector with 4 rows containing finger number
            mean_T1_D1(ind5,3) = ind5; % generates column vector with 4 rows containing layer number
            mean_T1_D1(ind5,4) = subject_code(ind);
            mean_T1_D1(ind5,5) = ind;
            mean_T1_D1(ind5,6) = vertices_d1_S1;
            mean_T1_D1(ind5,7) = h;

            mean_T1_D2(ind5,1) = nanmean(T1_D2(ind5,1:end));
            mean_T1_D2(ind5,2) = 2;
            mean_T1_D2(ind5,3) = ind5;
            mean_T1_D2(ind5,4) = subject_code(ind);
            mean_T1_D2(ind5,5) = ind;
            mean_T1_D2(ind5,6) = vertices_d2_S1;
            mean_T1_D2(ind5,7) = h;

            mean_T1_D3(ind5,1) = nanmean(T1_D3(ind5,1:end));
            mean_T1_D3(ind5,2) = 3;
            mean_T1_D3(ind5,3) = ind5;
            mean_T1_D3(ind5,4) = subject_code(ind);
            mean_T1_D3(ind5,5) = ind;
            mean_T1_D3(ind5,6) = vertices_d3_S1;
            mean_T1_D3(ind5,7) = h;

            mean_T1_all(ind5,1) = nanmean(T1_all(ind5,1:end));
            mean_T1_all(ind5,2) = 0;
            mean_T1_all(ind5,3) = ind5;
            mean_T1_all(ind5,4) = subject_code(ind);
            mean_T1_all(ind5,5) = ind;
            mean_T1_all(ind5,6) = vertices_d1_S1 + vertices_d2_S1 + vertices_d3_S1; %+vertices_d3_S1 + vertices_d4_S1 + vertices_d5_S1;
            mean_T1_all(ind5,7) = h;
         end
% 
         % write data to matfile
         % save(sprintf('%s_T1_mhand_DDlayer_Layers_mean_thrs_70.mat',subject{ind}),'mean_T1_D1');
         % save(sprintf('%s_T1_mface_DDlayer_Layers_mean_thrs_70.mat',subject{ind}),'mean_T1_D2');
         % save(sprintf('%s_T1_mfoot_DDlayer_Layers_mean_thrs_70.mat',subject{ind}),'mean_T1_D3');
         % save(sprintf('%s_T1_mface_mhan_mfoot_DDlayer_Layers_mean_thrs_70.mat',subject{ind}),'mean_T1_all');
% 
         table_allsubjects_fingers_long = [table_allsubjects_fingers_long; mean_T1_D1(:,:);mean_T1_D2(:,:);mean_T1_D3(:,:)];

         %l=length(table_allsubjects_fingers_short(ind,:))
         if h==1
            table_allsubjects_fingers_short(ind,:) = [transpose(mean_T1_D1(:,1)) transpose(mean_T1_D2(:,1)) transpose(mean_T1_D3(:,1)) subject_code(ind) ind nanmean(mean_T1_D1(:,6)) nanmean(mean_T1_D2(:,6)) nanmean(mean_T1_D3(:,6)) h];
         else
             table_allsubjects_fingers_short_na(ind,:) = [transpose(mean_T1_D1(:,1)) transpose(mean_T1_D2(:,1)) transpose(mean_T1_D3(:,1)) subject_code(ind) ind nanmean(mean_T1_D1(:,6)) nanmean(mean_T1_D2(:,6)) nanmean(mean_T1_D3(:,6)) h];
         end
% % 
         save(sprintf('%s_T1_mhand_DDlayer_mean_thrs_%s_%s_70.mat',subject{ind},hemi{h},res{r}),'mean_T1_D1');
         save(sprintf('%s_T1_mface_DDlayer_mean_thrs_%s_%s_70.mat',subject{ind},hemi{h},res{r}),'mean_T1_D2');
         save(sprintf('%s_T1_mfoot_DDlayer_mean_thrs_%s_%s_70.mat',subject{ind},hemi{h},res{r}),'mean_T1_D3');
         save(sprintf('%s_T1_mface_mhand_mfoot_DDlayer_mean_thrs_%s_%s_70.mat',subject{ind},hemi{h},res{r}),'mean_T1_all');

     else
        warning('Nothing done. Getting next subject...');
     end   

    end

    end

end
% C = [table_allsubjects_fingers_short table_allsubjects_fingers_short_na]
% C_21 = [table_allsubjects_fingers_short_21 table_allsubjects_fingers_short_21_na]
% 
% save(sprintf('allsubjects_T1_mFaceHandFoot_short_21Layers_%s_70.mat',res{1}),'C_21');
% save(sprintf('allsubjects_T1_mFaceHandFoot_long_21Layers_%s_70.mat',res{1}),'table_allsubjects_fingers_long_21');
% csvwrite(sprintf('allsubjects_T1_mFaceHandFoot_short_21L_%s_70mR.csv',res{1}),C_21);
% csvwrite(sprintf('allsubjects_T1_mFaceHandFoot_long_21L_%s_70mR.csv',res{1}),table_allsubjects_fingers_long_21);
% 
% % save('allsubjects_T1_mFaceHandFoot_short_DDLayers_70.mat','table_allsubjects_fingers_short');
% % save('allsubjects_T1_mFaceHandFoot_long_DDLayers_70.mat','table_allsubjects_fingers_long');
% % csvwrite('allsubjects_T1_mFaceHandFoot_short_DDLayers.csv',table_allsubjects_fingers_short);
% % csvwrite('allsubjects_T1_mFaceHandFoot_long_DDLayers.csv',table_allsubjects_fingers_long);
% 
% save(sprintf('allsubjects_T1_mFaceHandFoot_short_DDLayers_%s_70.mat',res{1}),'C');
% save(sprintf('allsubjects_T1_mFaceHandFoot_long_DDLayers_%s_70.mat',res{1}),'table_allsubjects_fingers_long');
% csvwrite(sprintf('allsubjects_T1_mFaceHandFoot_short_DDL_%s_70mR.csv',res{1}),C);
% csvwrite(sprintf('allsubjects_T1_mFaceHandFoot_long_DDL_%s_70mR.csv',res{1}),table_allsubjects_fingers_long);
