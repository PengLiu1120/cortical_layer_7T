function pl_layer_Gaussian_Tuning_Curve_Vertical(ModelFiles)
% Input order: Hemis, ModelFiles, ROIFiles

% Fits a 1D Gaussian tuning curve model with bars moving vertically

% .........................................................................
% Inputs
% .........................................................................
% Hemis: Hemispheres ('lh' and/or 'rh') [cell]
% ModelFiles: Cell array with SamSrf data files (without extension) [cell]
% ROIFiles: ROI label to restrict analysis [cell]

% Inputs are optional. If undefined, a dialog is opened for user selection.
%
% Wriiten by DSS
% Last updated 18 Feb 2022, by P.Liu

%% ........................................................................Standard 2D Gaussian pRF
% .........................................................................Set defaults
% .........................................................................Coarse fit only
Coarse_Fit_Only = [true false];
% .........................................................................Which pRF model function?
Model.Prf_Function = @(P,ApWidth) prf_gaussian_rf(0, P(1), P(2), ApWidth);
% .........................................................................File name to indicate type of pRF model
Model.Name = 'pTC_ver'; 
% .........................................................................Names of parameters to be fitted
% .........................................................................Mu: centre location, y, Sigma: size
Model.Param_Names = {'Mu_v'; 'Sigma'};
% .........................................................................Which of these parameters are scaled
% .........................................................................1: true, 0: false
Model.Scaled_Param = [1 1];
% .........................................................................Which parameters must be positive?
% .........................................................................1: true, 0: false
Model.Only_Positive = [0 1];
% .........................................................................Scaling factor of the stimulus space (e.g. eccentricity)
% .........................................................................No Scaling
Model.Scaling_Factor = 1;
% .........................................................................Temporal resolution of stimulus apertures (can be faster than scanner TR if downsampling predictions)
Model.TR = 2;
% .........................................................................HRF file or vector to use (empty = canonical) 
% .........................................................................If use individual HRF, HRF info should be converted to SamSrf format
Model.Hrf = [];
% .........................................................................Aperture file
Model.Aperture_File = 'aps_verbars_vec'; 

% Optional parameters
% .........................................................................Limit data to above certain noise ceiling?
% .........................................................................Correlation coefficient
Model.Noise_Ceiling_Threshold = 0; 
% .........................................................................If true, uses coarse fit for bad slow fits
Model.Replace_Bad_Fits = false; 
% .........................................................................If > 0, smoothes data for coarse fit
Model.Smoothed_Coarse_Fit = 0; 
% .........................................................................If true, only runs the coarse fit
% Model.Coarse_Fit_Only = false;
% .........................................................................Define a Srf file to use as seed map
Model.Seed_Fine_Fit = ''; 
% .........................................................................Define threshold for what to include in fine fit
Model.Fine_Fit_Threshold = 0.00; 
% .........................................................................Defines block size for coarse fit (reduce if using large search space)
Model.Coarse_Fit_Block_Size = 10000;
% .........................................................................If true, parameter 1 & 2 are polar (in degrees) & eccentricity coordinates
Model.Polar_Search_Space = false; 
% .........................................................................Use for microtime resolution if stimulus timing is faster than TR
% Model.Downsample_Predictions = 10;

% Search grid for coarse fit
% .........................................................................Mu search grid
Model.Param1 = -1.05 : 0.05 : 1.05;
% .........................................................................Sigma search grid
Model.Param2 = 2 .^ (-5.6 : 0.1 : 1);
% .........................................................................
Model.Param3 = 0; % Unused
Model.Param4 = 0; % Unused
Model.Param5 = 0; % Unused

%% ........................................................................Fitting process

for i_fit = 1:size(Coarse_Fit_Only, 2)
    
    Model.Coarse_Fit_Only = Coarse_Fit_Only(i_fit);
    
    CurrSrfFile = ['vol_' ModelFiles];
    
    MapFile = samsrf_fit_prf(Model, CurrSrfFile);
    
end

end