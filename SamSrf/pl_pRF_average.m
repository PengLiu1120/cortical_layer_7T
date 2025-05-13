function pl_pRF_average(Hemis,Run)
% Input order: Hemis, Run

% Average files first converted by mgh2srf from FreeSurfer MGH format then
% concatenated by pl_pRF_concatenate

% .........................................................................
% Inputs
% .........................................................................
% Hemis: Hemispheres ('lh' and/or 'rh') [cell]
% Run: The run number ('01' and/or '02') [cell]
% .........................................................................
% Outputs
%..........................................................................
% Note: Noise ceiling can only be performed with even run numbers
% .........................................................................

% Written by P.Liu
% Last Updated 20 Jan 2022 by P.Liu
%% ........................................................................Function

% .........................................................................Loop through hemi
for i_hemi = 1:size(Hemis, 2)
    
    CurrHemi = Hemis{i_hemi};
    
    Functional = [];
    CurrFunctional = [];
    CurrData = [];
    Data = [];
    
    % .....................................................................Loop through runs
    for i_run = 1:size(Run,2)
        
        CurrRun = Run{i_run};
        
        CurrFile_info = dir([CurrHemi '_' CurrRun '.mat']);
        CurrFile = load(CurrFile_info.name);
        
        CurrFunctional = CurrFile_info.folder;
        Functional = [Functional; {fullfile(CurrFunctional, [CurrHemi '_' CurrRun])}];
        CurrData{i_run} = CurrFile.Srf.Data;
        
        Data = cat(3, Data, CurrData{i_run});
        
    end
    
    Srf = CurrFile.Srf;
    Srf.Functional = Functional;

    Srf.Data = Data;
    
    % .....................................................................Calculate Noise Ceiling
    OddRuns = nanmean(Srf.Data(:,:,1:2:end), 3);
    EvenRuns = nanmean(Srf.Data(:,:,2:2:end), 3);
    
    % .....................................................................Loopthrough vertices
    Srf.Noise_Ceiling = NaN(1, size(Srf.Data,2));
    
    for v = 1:size(Srf.Data, 2)
        
        % .................................................................Correlation between odd & even runs
        Rho_xxp = corr(OddRuns(:,v), EvenRuns(:,v));
        % .................................................................Spearman-Brown prediction formula
        Srf.Noise_Ceiling(v) = (2*Rho_xxp) / (1+Rho_xxp);
        % .................................................................For observable correlation we would need to take square root so without the square root this is R^2!
        
    end
    
    % .....................................................................Calculate mean across runs
    Srf.Data = nanmean(Srf.Data, 3);
    
    FileName = [CurrHemi '_Average'];
    save(FileName, 'Srf', '-v7.3');
    
end

end