function pl_layer_pRF_concatenate(Order, Run)
% Input order: Order, Run

% Concatenate files converted by mgh2srf from FreeSurfer MGH format

% .........................................................................
% Inputs
% .........................................................................
% Order: The stimulation order ('forward' and/or 'backward') [cell]
% Run: The run number ('01' and/or '02') [cell]
% .........................................................................
% Outputs
%..........................................................................

% Written by P.Liu
% Last Updated 20 Jan 2022
%% ........................................................................Function    
% .........................................................................Loop through runs
for i_run = 1:size(Run,2)
    
    CurrRun = Run{i_run};
    
    Functional = [];
    Data = [];
    FileName = [];
    
    % .....................................................................Loop through order
    for i_order = 1:size(Order, 2)
        
        CurrOrder = Order{i_order};
        CurrFile = load(['vol_' CurrOrder CurrRun]);
        
        Functional = [Functional; CurrFile.Srf.Functional];
        Data = [Data; CurrFile.Srf.Data];
        
        FileName = [FileName CurrOrder CurrRun];
        
    end
    
    Srf = CurrFile.Srf;
    Srf.Functional = Functional;
    Srf.Data = Data;
    
    Srf.Values = {};
    
    for i_vol = 1:size(Srf.Data,1)
        
        Srf.Values{i_vol,1} = ['Volume #' num2str(i_vol)];
        
    end
    
    FileName = ['vol_' FileName];
    save(FileName, 'Srf', '-v7.3');
    
end

end