% Step 1: Calculating Resting-State Fluctuation Amplitude (RSFA)
% RSFA = The temporal SDs of the corresponding timeseries

%% ........................................................................Set Defaults
FileName = [];

%% ........................................................................RSFA Calculation

resting = spm_vol('sub-ajz_resting_state_registered_to_ajz_run-01_T1map.nii');

resting_state = spm_read_vols(resting);

RSFA = std(resting_state);

%% ........................................................................BOLD Normalisation

for i_run = 1:size(Run,2)
    
    CurrRun = Run{i_run};
    
    for i_order = 1:size(Order, 2)
        
        CurrOrder = Order{i_order};
        
        CurrFile = load(['vol_' CurrOrder CurrRun]);
        
        FileName = [FileName CurrOrder CurrRun];
        
        Data = CurrFile.Srf.Data;
        
        NormData = zeros(size(Data));
        
        for i_row = 1:size(Data,1)
            
            for i_col = 1:size(Data,2)
                
                NormData(i_row,i_col) = Data(i_row,i_col)-RSFA(i_col)*(1-Data(i_row,i_col)^((0.2-1)/0.38));
                
            end
            
        end
        
    end
    
end