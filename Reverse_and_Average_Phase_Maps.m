%function AvgPhaseMaps
%
% Select the two time series you want to reverse and average. 
% You can select more than one file in each dialog,
% but be sure you separate them by their direction!
%

% definitions
datadir = '/Users/pliu/Documents/DataAnalysis0512/Dataset/lym720/0624/D3/down';
% script_dir = '';
x_dim = 288;
y_dim = 288;
no_slices = 26;
volumes = 256;
% names = {'BI3T/Fourier/' 'DJ6T/Fourier/' 'EC7T/Fourier/' 'FMFT/Fourier/' 'GD5T/Fourier/' 'HJJT/Fourier/' 'KSWT/Fourier/' 'NL5T/Fourier/' 'NW3T/Fourier/' 'OK2T/Fourier/' 'RSIT/Fourier/' 'SC1T/Fourier/' 'SLJT/Fourier/' 'ZK4T/Fourier/'};

% specifies spm path (for mac)
% addpath /Applications/spm8

% some definitions
m.mean = [];
n = 1;

% start of loop (one for each subject)
% for ind = 1:length(names)

% go to subject folder
cd (datadir)
% subjectdir = sprintf('%s',names{ind});
% cd (subjectdir)

% define files
% fname_D2_D5 = sprintf('%sD2_D5_rw_audata_cut.nii',names{ind});
% fname_D5_D2 = sprintf('%sD5_D2_rw_audata_cut.nii',names{ind});

% Data2Read_D2_D5 = fullfile(datadir,fname_D2_D5); %creates string character
% Data2Read_D5_D2 = fullfile(datadir,fname_D5_D2);
% 
% % read header (needs SPM toolbox)
% HeaderInfo_D2_D5 = spm_vol(Data2Read_D2_D5); % creates struc file with header info, including .dim
% HeaderInfo_D5_D2 = spm_vol(Data2Read_D5_D2);
% 
% % use spm_read_vols to read in the data
% % size(Data)=[128,128,44,256]
% Data_D2_D5 = spm_read_vols(HeaderInfo_D2_D5);
% Data_D5_D2 = spm_read_vols(HeaderInfo_D5_D2);

%% Reverse and save anticlockwise scan (name: D5_D2_rw_audata_cut_reversed.nii)
% Data_back_D5_D2 = flip(Data_D5_D2,4);
% 
% for i=1:volumes
%     HeaderInfo_D5_D2(i).fname = ['D5_D2_rw_audata_cut_reversed.nii'];
%     spm_write_vol(HeaderInfo_D5_D2(i), Data_back_D5_D2(:,:,:,i));
% end
% 
% % extract values from second and third volumes of reversed scan, 
% % and average second and third volume
% Volume2 = NaN(x_dim,y_dim,no_slices);
% Volume3 = NaN(x_dim,y_dim,no_slices);
% Volume2 = Data_back_D5_D2(1:x_dim,1:y_dim,1:no_slices,2);
% Volume3 = Data_back_D5_D2(1:x_dim,1:y_dim,1:no_slices,3);
% Volume_avg = (Volume2+Volume3)/2;
% 
% % replace first volume by average, 
% % and delete last volume of reversed scan (name: D5_D2_rw_audata_cut_reversed_del.nii)
% Data_back_D5_D2_del = NaN(x_dim,y_dim,no_slices);
% Data_back_D5_D2_del(1:x_dim,1:y_dim,1:no_slices,1) = Volume_avg;
% Data_back_D5_D2_del(1:x_dim,1:y_dim,1:no_slices,2:volumes) = Data_back_D5_D2(1:x_dim,1:y_dim,1:no_slices,1:volumes-1);
% 
% for i=1:volumes
%     HeaderInfo_D5_D2(i).fname = ['D5_D2_rw_audata_cut_reversed_del.nii'];
%     spm_write_vol(HeaderInfo_D5_D2(i), Data_back_D5_D2_del(:,:,:,i));
% end


%% Average D2_D5 and D5_D2 images, and save (name: avg_rw_audata.nii)
% for j=1:volumes
%     Volume1 = NaN(x_dim,y_dim,no_slices);
%     Volume2 = NaN(x_dim,y_dim,no_slices);
%     Data_avg = NaN(x_dim,y_dim,no_slices,volumes);
%     
%     Volume1 = Data_D2_D5(1:x_dim,1:y_dim,1:no_slices,j);
%     Volume2 = Data_back_D5_D2_del(1:x_dim,1:y_dim,1:no_slices,j);
%     Volume_avg = (Volume1+Volume2)/2;
%     
%     Data_avg(1:x_dim,1:y_dim,1:no_slices,j)=Volume_avg;
%     HeaderInfo_D2_D5(j).fname = ['avg_rw_audata.nii'];
%     spm_write_vol(HeaderInfo_D2_D5(j), Data_avg(:,:,:,j));
% end


% n = n + 1
    
% end

%cd (script_dir)
%% 
fname_down = sprintf('adata2.nii');
Data2Read_down = fullfile(datadir,fname_down);
HeaderInfo_down = spm_vol(Data2Read_down);
Data_down = spm_read_vols(HeaderInfo_down);

Data_back_down = flip(Data_down,4);

for i=1:volumes
    HeaderInfo_down(i).fname = ['adata2.nii'];
    spm_write_vol(HeaderInfo_down(i), Data_back_down(:,:,:,i));
end

% extract values from second and third volumes of reversed scan, 
% and average second and third volume
Volume2 = NaN(x_dim,y_dim,no_slices);
Volume3 = NaN(x_dim,y_dim,no_slices);
Volume2 = Data_back_down(1:x_dim,1:y_dim,1:no_slices,2);
Volume3 = Data_back_down(1:x_dim,1:y_dim,1:no_slices,3);
Volume_avg = (Volume2+Volume3)/2;

% replace first volume by average, 
% and delete last volume of reversed scan
Data_back_down_reversed = NaN(x_dim,y_dim,no_slices);
Data_back_down_reversed(1:x_dim,1:y_dim,1:no_slices,1) = Volume_avg;
Data_back_down_reversed(1:x_dim,1:y_dim,1:no_slices,2:volumes) = Data_back_down(1:x_dim,1:y_dim,1:no_slices,1:volumes-1);

for i=1:volumes
    HeaderInfo_down(i).fname = ['adata2_reversed.nii'];
    spm_write_vol(HeaderInfo_down(i), Data_back_down_reversed(:,:,:,i));
end