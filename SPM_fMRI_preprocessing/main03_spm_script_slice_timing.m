%__________________________________________________________________________
% Batch script:
%   - this script performs slice timing of an fMRI time series
% modified August 2018 by Esther Kuehn (esther.kuehn@dzne.de)
%
%__________________________________________________________________________


%===========================================================================
% user specified settings
%===========================================================================
clear all

%===========================================================================
% user specified path to the SPM8 directory
% please check SPM and MATLAB version
%===========================================================================
%
% To set the path, use the command 'SPM' before you start 'matlab' 

n_sess          = 1;      % number of sessions
n_scans         = 256;    % scans per session (294 for real, 550 for observed)
TR              = 2.0;      % repetition time in s
num_slices      = 36;     % number of slices
ref_slice       = 18;     % reference slice for slice timing (e.g. first slice or middle slice)
			  % IMPORTANT: The microtime-onset needs to be adjusted according to the reference slice!
			  %            The default is SPM.xBF.T0 = 8 (middle slice), but if you choose the first slice
			  %	       as reference you have to change SPM.xBF.T0 = 1. 
			  %	       The parameter can be find in the parameter-block for defining the design matrix.
              
dir_base        = '/media/...'; % directory which contains all subject-specific directory
dir_functional  = ''; % directory of functional data (NIfTI files)
name_scans      = 'data';       % name of the NIfTI functional data-file 
sess_prfx       = 'sess';       % prefix of session directories (only needed for multi-session analyses)

delete_files	= 0;  % delete intermediate steps during pre-processing? 1=yes, 0=no
data_scale      = 1; % scale data before realign and unwarp: 
                      % data scaling factor, ENSURE THAN MAXUMUM VALUE WILL BE < 32000
slice_timing    = 2;  % slice timing? 2=before realignment, 1=after realignment, 0=no correction  (because interleaved 			acquisition)
slice_order     = 0;  % slice order? 2=custom, 1=ascending, 0=descending
slice_ord       = [num_slices:-1:1];	% Syntax: [first_slice:increment/decrement:last_slice]
							% Examples:
							% descending acquisition [num_slices:-1:1]
							% ascending acquisition [1:1:num_slices]
							% interleaved acquisition (ascending - even slices first): 									[2:2:num_slices 1:2:num_slices-1]
  
start_analysis  = 0;  % start of analysis (0 means run normally)

% names of subjects (list all names of subject directories), example
% subject, you can add more subjects
name_subj = {'ulb827/D1_D5'};



%===========================================================================
% end of user specified settings
%===========================================================================



% call SPM to have the graphics windows in place
%===========================================================================
spm fmri
global defaults
defaults.stats.maxmem   = 2^35;
fs = filesep; % platform-specific file separator
k  = 1;

% set slicetime parameters
if slice_timing>0
      
    slice_time = TR /num_slices;
    slice_gap  = slice_time;	    % continuous scanning
    
    % descending slice acquisition
    if slice_order == 0
         slice_ord = [num_slices:-1:1];
    end 

    % ascending slice acquisition
    if slice_order == 1
         slice_ord = [1:1:num_slices];
    end

end    


% loop through all subjects
%===========================================================================
while (k<=length(name_subj)),
    
    fprintf(1,'==================================\n');
    fprintf(1,'Starting analysis for subject %s\n',name_subj{k});
    fprintf(1,'==================================\n');

        
    % define subject-specific scan directories
    % ----------------------------------------
    if n_sess > 1 
    	for sess    = 1:n_sess,
    	    dir_scans{sess} = [dir_base fs name_subj{k} fs dir_functional fs sess_prfx num2str(sess)];
    	end
    else
    	dir_scans{1} = [dir_base fs name_subj{k} fs dir_functional];
    end
    
    % slice timing
    %===========================================================================
    if start_analysis <= 1
        % Go through all sessions of the subject
    	%---------------------------------------
    	for sess = 1:n_sess,
	
	
            % data scaling as first step
	    %---------------------------
	    if data_scale > 1.1
	       cd (dir_scans{sess})
	       copyfile([name_scans '.nii'],[name_scans '_unscaled.nii']);
	       if (data_scale > 50)
	          warning('Check size of scaling factor');
	       end
	       
	       for nsc = 1:n_scans
	          P   = fullfile(dir_scans{sess},[name_scans '.nii,' mat2str(nsc)]); 
		  Vsc = spm_vol(P);
		  X   = spm_read_vols(Vsc); 
		  X   = X .* data_scale;
		  Vsc = spm_write_vol(Vsc,X); 
	       end	    
	    end
	
	
    	    % Slice-timing correction 
    	    %------------------------
	    prefix = '';
    	    if slice_timing == 2
	    
	        % load NIfTI files
    		%-----------------------------
		prefix = 'a';
    		cd (dir_scans{sess})
		P = fullfile(dir_scans{sess},[name_scans '.nii']);
     	 
    		% Slicetime correction
    		%------------------------------------
    		spm_slice_timing(P, slice_ord, ref_slice, [slice_time slice_gap], 'a');

    	    end
    	    end
    	    end
       
    
    
    
    % Switch to next subject
    %=======================================
    k   = k + 1;
    clear SPM P tmpP;

end     % end of main loop

cd (dir_base);
