function ss_samsrf_mgh2srf(Hemis, varargin)
% Input order: Hemis, Combi, FldSrf, FuncPattern, Norm, Avg, NoiseCeiling,
% AnatPath
%
% Converts FreeSurfer MGH file to SamSrf surface format.
%
% ----------------------------------------------------------------------------
% Inputs
% ----------------------------------------------------------------------------
% Hemis             - Hemispheres('lh' and/or 'rh') [cell]
% Combi             - Type of data combination ('runwise' or 'accrun') [char]
% FldSrf            - 'Surf' folder name [char]
% FuncPattern       - Pattern of functional file name [char]
% Norm              - Toggles whether normalization should be applied or not
%                     [logical]
% Avg               - Toggles whether multiple runs will be averaged (true) or
%                     concatenated (false) [logical]
% NoiseCeiling      - Noise ceiling (true=calculates noise ceiling) [logical]
% AnatPath          - Anatomy path [char]
% ----------------------------------------------------------------------------
% Outputs
% ----------------------------------------------------------------------------
% -/-
% ----------------------------------------------------------------------------
% 27/12/2020: Generated (DSS, SS)
% 21/01/2022: Last modified (PL)

%% ...........................................................................Some defaults

DefInputs = {'runwise' '..\..\surf\' 'ubf*.nii'  true false 'false' '../anatomy/'};
% Assign user-specifie inputs
DefInputs(1:length(varargin)) = varargin;
[Combi, FldSrf, FuncPattern, Norm, Avg, NoiseCeiling, AnatPath] = DefInputs{:};

%% ...........................................................................List functional files

FuncList = dir(FuncPattern);

%% ...........................................................................Check existence of functional files

if ~isempty(FuncList)
    
    %% .......................................................................Get names of functional files and transpose
    
    FuncNames = {FuncList.name}';
    
    %% .......................................................................Loop through functional files
    
    for i_fu_n = 1:length(FuncNames)
        
        %% ...................................................................Cut functional file name
       
        FuncNames{i_fu_n} = FuncNames{i_fu_n}(1:end-4);
        
    end
    
    %% .......................................................................Loop through hemis
    
    for i_hemi = 1:size(Hemis, 2)
        
        %% ...................................................................Which data combination mode?
        
        CurrHemi = Hemis{i_hemi};
        CurrHemiSurf = fullfile(FldSrf,  CurrHemi);
        CurrFuncNames = FuncNames(contains(FuncNames, CurrHemi));
        
        switch Combi
            case 'runwise' % Run-wise projection 
                
                %% ...........................................................Loop through functional files
                
                for i_fu_l = 1:length(CurrFuncNames)
                    samsrf_mgh2srf(fullfile(FuncList(i_fu_l).folder, CurrFuncNames{i_fu_l}), ...
                        CurrHemiSurf, Norm, Avg, NoiseCeiling, AnatPath)
                end
            
            case 'accrun' % Across run projection  
                samsrf_mgh2srf(fullfile(FuncList(1).folder, CurrFuncNames), ...
                    CurrHemiSurf, Norm, Avg, NoiseCeiling, AnatPath)
        end
    end
    
else
    error('No files found. Try again.')
end

end