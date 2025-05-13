function pl_copylabel(SubjPath, PRFPath, Subject, Label, Loc)
% function pl_copylabelregister(Hemis, SubjPath, SessPath, SPMPath, PRFPath, Subject, Label, Condition, Image, Fourier, FuncSession)
% Input order: Hemis, ModelFiles, PathSPM, PathPRF

% Copy labels used for data extraction
% Labels are 3b hand + localiser overlap

% Only left hemisphere

% .........................................................................
% Inputs
% .........................................................................
% Hemis: Hemispheres ('lh' and/or 'rh') [cell]
% SubjPath: Freesurfer subject path [cell]
% SessPath: Freesurfer session path [cell]
% SPMPath: Path to spm_ folder [cell]
% PRFPath: Path to prf_ folder [cell]
% Subject: Subject [cell]
% Label: label [cell]
% Condition: The stimulation condition (D2 and/or D3 and/or D2+D3) [cell]
% Image: 'image' [cell]
% Fourier: '_Fourier' [cell]
% FuncSession: '_backward_1' [cell]
% .........................................................................
% Outputs
% .........................................................................
% Written by P.Liu
% Last updated 25 Mar 2022 by P.Liu
%% ........................................................................Function
% .........................................................................copy label

LabelDir = fullfile(SubjPath,Subject,Label);

% Hand_Loc = cell2mat(['lh.3b_hand_Loc_' Condition '.label']);
% BA3b_Loc = cell2mat(['lh.3b_Loc_' Condition '.label']);
Localiser = ['lh.' Loc '_Loc.label'];

% Hand_New = cell2mat(['lh.hand_loc_' Condition '.label']);

% copyfile(strcat(LabelDir,'/',Hand_Loc), strcat(PRFPath,'/',Hand_New))
copyfile(fullfile(LabelDir,Localiser), PRFPath)

end