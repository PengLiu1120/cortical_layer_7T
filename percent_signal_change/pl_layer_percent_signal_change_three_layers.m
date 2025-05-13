%% Layer-specific Percent Signal Response
% .........................................................................
% This script calculate layer-specific percent signal change for functional
% MRI data.
% .........................................................................
% To estimate map amplitudes (in percentage), we started with the discrete
% Fourier transform response amplitude (hypotenuse given real and imaginary
% values) for each vertex within our ROIs. This value was multiplied by 2
% to account for positive and negative frequencies, again multiplied by 2
% to estimate peak-to-peak values, divided by the number of time points
% over which averaging was performed (to normalize the discrete FT
% amplitude), and divided by the average brightness of the functional
% dataset (excluding air). Finally, the value was multiplied by 100 for
% percentage response amplitude.
% Kuehn et al. 2018
% .........................................................................
% Functional data should first use Csurf to calculate Fourier map
% Take results file end with _x and mask with chosen mask, write out ROI
% with label.
% .........................................................................
% Concatenate all runs (i.e., forward runs and backward runs), using
% fslmaths to calculate mean image.
% Use fsl to skull strip air from the mean image, threshold 0.1.
% .........................................................................
% Inputs
% .........................................................................
% _x and _y values
% skull stripped mean image
% timepoints (i.e., number of TRs)
% .........................................................................
%                 hypot(_x,_y)
% %response = -------------------- * 2 * 2 * 100
%              timpoints * rawavg
% .........................................................................
% Outputs
% .........................................................................
% % signal response
% .........................................................................
% Written by P.Liu
% Optimized by P.Liu
% Email: peng.liu@med.ovgu.de
% Last modified by P.Liu 25 Apr 2022
%% ........................................................................Tidy up
clear all
close all
clc

%% ........................................................................Specify parameters
% .........................................................................Specify root path
RootPath = '/media/pliu/LayerPRF/LayerMapping/08_LayerExtraction/Fourier';

% .........................................................................Specify subjects
Subjects = {'ajz367' 'bkn792' 'bmg520' 'cxc075' 'czg996' 'frj712' 'ggp057' 'gph998' 'gxo876' 'hby152' 'ijt563' 'iwq192' 'kdy341' 'llh150' 'lpr469' 'nhm378' 'oms448' 'qet940' 'qxo538' 'sst050' 'unk742'};
Subj = {'ajz' 'bkn' 'bmg' 'cxc' 'czg' 'frj' 'ggp' 'gph' 'gxo' 'hby' 'ijt' 'iwq' 'kdy' 'llh' 'lpr' 'nhm' 'oms' 'qet' 'qxo' 'sst' 'unk'};

% .........................................................................Specify condtions, D2 and D2+D3 or D3 and D2+D3
Conditions = {'D2+D3' 'D2' 'D3'};
% 'D2+D3' 'D2' 'D3'

% .........................................................................Specify layers, deep and superficial
Layers = {'Deep' 'Middle' 'Superficial'};
layers = {'deep' 'middle' 'superficial'};

% .........................................................................Parameters
SignalChange = '%signalchange';
timepoints = 512; % Number of TRs

%% ........................................................................Calculate % signal resposne
for i_cond = 1:size(Conditions,2)
    
    CurrCond = Conditions{i_cond};
    
    for i_sub=1:size(Subjects, 2)
        
        CurrSubj = Subjects{i_sub};
        CurrS = Subj{i_sub};
        
        % .....................................................................Direct to subject path
        CurrSubjPath = fullfile(RootPath,CurrCond,CurrSubj);
        
        cd(CurrSubjPath)
        
        for i_layer = 1:size(Layers,2)
            
            CurrLayer = Layers{i_layer};
            Currlayer = layers{i_layer};
            
            % .............................................................Calculate average brightness
            % BrightnessImg = ['sub-' CurrS '_' CurrCond '_afourier_mean_registered_to_' CurrS '_run-01_T1map_brain.nii.gz'];
            % gunzip(BrightnessImg);
            Brightness_mean = spm_vol(['sub-' CurrS '_' CurrCond '_afourier_mean_registered_to_' CurrS '_run-01_T1map_brain.nii']);
            Fourier_mean = spm_read_vols(Brightness_mean);
            
            % .............................................................Replace 0s with NaN
            Fourier_mean(Fourier_mean == 0) = NaN;
            
            % .............................................................Calculate average brightness for each slice
            Fourier = [];
            for i = 1:size(Fourier_mean,3)
                Fourier_slice = Fourier_mean(:,:,i);
                Fourier_slice_avg = nanmean(Fourier_slice,'all');
                Fourier = [Fourier; Fourier_slice_avg];
            end
            
            % .............................................................Calculate average brightness across slice
            rawavg = nanmean(Fourier);
            
            % .............................................................Read out _x and _y values for overlap vertices
            x_file = load([CurrSubj '_Fourier_i_' CurrLayer]);
            x_name = [Currlayer '_layers'];
            x = x_file.(x_name);
            
            x(x == 0) = NaN;
            
            y_file = load([CurrSubj '_Fourier_r_' CurrLayer]);
            y_name = [Currlayer '_layers'];
            y = y_file.(y_name);
            
            y(y == 0) = NaN;
            
            percent_signal_change = [];
            
            Hyp = hypot(x,y);
            Hyp_mean = nanmean(Hyp);
            
            % .............................................................Calculate percent signal change
            percent_response = (Hyp_mean*2*2*100)/(timepoints*rawavg);
            
            percent_signal_change = [percent_signal_change; percent_response];
            
            FileName = [CurrCond '_' CurrLayer '_percentsignalchange'];
            save(FileName, 'percent_signal_change');
            
        end
        
    end
    
end