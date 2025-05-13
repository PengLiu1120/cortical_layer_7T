function pl_layer_pRF_average(Run)
% Input order: Run

% Average files first converted by mgh2srf from FreeSurfer MGH format then
% concatenated by pl_pRF_concatenate

% .........................................................................
% Inputs
% .........................................................................
% Run: The run number ('01' and/or '02') [cell]
% .........................................................................
% Outputs
%..........................................................................
% Note: Noise ceiling can only be performed with even run numbers
% .........................................................................

% Written by P.Liu
% Optimized by P.Liu
% Email: peng.liu@med.ovgu.de
% Last Updated 11 Oct 2022 by P.Liu
%% ........................................................................Function
Functional = [];
CurrFunctional = [];
CurrData = [];
Data = [];

% .........................................................................Loop through runs
for i_run = 1:size(Run,2)
    
    CurrRun = Run{i_run};
    
    CurrFile_info = dir(['vol_' CurrRun '.mat']);
    CurrFile = load(CurrFile_info.name);
    
    CurrFunctional = CurrFile_info.folder;
    Functional = [Functional; {fullfile(CurrFunctional, ['vol_' CurrRun])}];
    CurrData{i_run} = CurrFile.Srf.Data;
    
    Data = cat(3, Data, CurrData{i_run});
    
end

Srf = CurrFile.Srf;
Srf.Functional = Functional;

Srf.Data = Data;

% .........................................................................Calculate Noise Ceiling
OddRuns = nanmean(Srf.Data(:,:,1:2:end), 3);
EvenRuns = nanmean(Srf.Data(:,:,2:2:end), 3);

% .........................................................................Loopthrough vertices
Srf.Noise_Ceiling = NaN(1, size(Srf.Data,2));

for v = 1:size(Srf.Data, 2)
    
    % .....................................................................Correlation between odd & even runs
    Rho_xxp = corr(OddRuns(:,v), EvenRuns(:,v));
    % .....................................................................Spearman-Brown prediction formula
    Srf.Noise_Ceiling(v) = (2*Rho_xxp) / (1+Rho_xxp);
    % .....................................................................For observable correlation we would need to take square root so without the square root this is R^2!
    
end

% .........................................................................Calculate mean across runs
Srf.Data = nanmean(Srf.Data, 3);

FileName = 'vol_Average';
save(FileName, 'Srf', '-v7.3');

end