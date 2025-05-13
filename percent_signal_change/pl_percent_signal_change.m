%% % Signal Response
% .........................................................................
% This script calculate % signal change for functional MRI data.
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
RootPath = '/Volumes/IKND/AG_Kuehn/Peng/LayerPRF/LayerPRF';

% .........................................................................Specify subjects
Subjects = {'ajz367' 'bkn792' 'bmg520' 'cxc075' 'czg996' 'frj712' 'ggp057' 'gph998' 'gxo876' 'hby152' 'ijt563' 'iwq192' 'kdy341' 'llh150' 'lpr469' 'nhm378' 'oms448' 'qet940' 'qxo538' 'sst050' 'unk742'};
% 'ajz367' 'bkn792' 'bmg520' 'cxc075' 'czg996' 'frj712' 'ggp057' 'gph998' 'gxo876' 'hby152' 'ijt563' 'iwq192' 'kdy341' 'llh150' 'lpr469' 'nhm378' 'oms448' 'qet940' 'qxo538' 'sst050' 'unk742'
Subj = {'ajz' 'bkn' 'bmg' 'cxc' 'czg' 'frj' 'ggp' 'gph' 'gxo' 'hby' 'ijt' 'iwq' 'kdy' 'llh' 'lpr' 'nhm' 'oms' 'qet' 'qxo' 'sst' 'unk'};
% 'ajz' 'bkn' 'bmg' 'cxc' 'czg' 'frj' 'ggp' 'gph' 'gxo' 'hby' 'ijt' 'iwq' 'kdy' 'llh' 'lpr' 'nhm' 'oms' 'qet' 'qxo' 'sst' 'unk'

% .........................................................................Specify condtions, D2 and D2+D3 or D3 and D2+D3
Conditions = {'D2' 'D3' 'D2+D3'};

% .........................................................................Specify localisers, hand and 3b
Loc = {'3b_hand'};

% .........................................................................Parameters
Results = 'Results';
SignalChange = '%signalchange';
timepoints = 512; % Number of TRs

%% ........................................................................Calculate % signal resposne
for i_sub=1:size(Subjects, 2)
    
    CurrSubj = Subjects{i_sub};
    CurrS = Subj{i_sub};
    
    % .....................................................................Direct to subject path
    CurrSubjPath = fullfile(RootPath,CurrSubj);
    
    for i_cond = 1:size(Conditions,2)
        
        CurrCond = Conditions{i_cond};
        
        % .................................................................Direct to % signal change path
        CurrResultPath = fullfile(CurrSubjPath,Results,SignalChange,CurrCond);
        
        cd(CurrResultPath)
        
        % .................................................................Calculate average brightness
        BrightnessImg = 'afourier_mean_brain.nii.gz';
        gunzip(BrightnessImg);
        Brightness_mean = spm_vol('afourier_mean_brain.nii');
        Fourier_mean = spm_read_vols(Brightness_mean);
        
        % .................................................................Replace 0s with NaN
        Fourier_mean(Fourier_mean == 0) = NaN;
        
        % .................................................................Calculate average brightness for each slice
        Fourier = [];
        for i = 1:size(Fourier_mean,3)
            Fourier_slice = Fourier_mean(:,:,i);
            Fourier_slice_avg = nanmean(Fourier_slice,'all');
            Fourier = [Fourier; Fourier_slice_avg];
        end
        
        % .................................................................Calculate average brightness across slice
        rawavg = nanmean(Fourier);

        for i_loc = 1:size(Loc,2)
            
            CurrLoc = Loc{i_loc};
            
            % .............................................................Specify ROI label
            Label = ['lh.' CurrLoc '_' CurrCond '.label'];
            
            % .............................................................Read out _x and _y values for overlap vertices
            vertex = load(Label);
            format long;
            Vertices = vertex(:,1);
            Coordinates_X = vertex(:,2);
            Coordinates_Y = vertex(:,3);
            Coordinates_Z = vertex(:,4);
            X = vertex(:,5);
            Y = vertex(:,6);
            
            percent_signal_change = [];
            
            % .............................................................Calculate hypotense for each vertex
            for i_ver = 1:size(X)
                
                Hyp = hypot(X(i_ver),Y(i_ver));
                
                % .........................................................Calculate percent signal change
                percent_response = (Hyp*2*2*100)/(timepoints*rawavg);
                
                percent_signal_change = [percent_signal_change; percent_response];
                
                extra_col = zeros(size(X));
                
            end
            
            Hyp_avg = hypot(X,Y);
            Hyp_avg_mean = nanmean(Hyp_avg);
            percent_signal_change_avg = (Hyp_avg_mean*2*2*100)/(timepoints*rawavg);
            
            FileName = [CurrLoc '_' CurrCond '_percentsignalchange'];
            save(FileName,'Vertices','percent_signal_change','percent_signal_change_avg');
            
            response = [Vertices Coordinates_X Coordinates_Y Coordinates_Z percent_signal_change extra_col];
            responseName = ['lh.' CurrLoc '_' CurrCond '_percentsignalchange.label'];
            dlmwrite(responseName, response, 'delimiter','\t');
            
        end
        
    end
    
end