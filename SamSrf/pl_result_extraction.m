function pl_result_extraction(Condition, Loc)
% function pl_result_extraction(Hemis,Condition)
% Input Order: Hemis, Conditions

% Extract results based on regions of interest from modelling fit

% .........................................................................
% Inputs
% .........................................................................
% Hemis: Hemispheres ('lh' and/or 'rh') [cell]
% Conditions: Stimulation conditions (D2, D3 and D2+D3) [cell]
% .........................................................................
% Outputs
%..........................................................................

% Written by P.Liu
% Last Updated 11 Apr 2022 by P.Liu
%% ........................................................................Function

% .........................................................................Loop through Hemis
% for i_hemi = 1:size(Hemis, 2)

% CurrHemi = Hemis{i_hemi};

% .........................................................................Load result .mat file
% load([CurrHemi '_pTC_ver.mat']);
load('lh_pTC_ver.mat');

% .........................................................................Expand Srf
Srf = samsrf_expand_srf(Srf);
Values = Srf.Values;

% .........................................................................Load ROI label file
% Hand_Loc = samsrf_loadlabel([CurrHemi '.3b_hand_loc_' Condition]);
ROI_Loc = samsrf_loadlabel(['lh.' Loc '_Loc']);

% .........................................................................Mask label onto data
Data_ROI = Srf.Data(:,ROI_Loc);

indexR2 = strcmp(Srf.Values, 'R^2');

TrueData = Data_ROI;

TrueData(:,TrueData(indexR2,:) == 0) = [];

% FileName = [CurrHemi '_' Condition '_model'];
FileName_Hand = ['lh_' Loc '_' Condition '_model'];
save(FileName_Hand, 'TrueData', 'Values', '-v7.3');
    
% end

end