%__________________________________________________________________________
% Batch script:
%   - preprocessing
%   - GLM estimation
%   - contrast definition & evaluation
% for a series of analyses from different subjects.
% Modified by Burkhard, Karsten & Joeran MPI, January 2010
%
% August 2010: Added pre-orient all images with the MNI space by co-registration with T1 template, KM
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
% To set the path, use the command 'SPM' before you start 'matlab7' 

n_sess          = 1;      % number of sessions
n_scans         = 264;    % scans per session
TR              = 2.0;    % repetition time in s
num_slices      = 30;     % number of slices
ref_slice       = 2;      % reference slice for slice timing (e.g. first slice or middle slice)
			  % IMPORTANT: The microtime-onset needs to be adjusted according to the reference slice!
			  %            The default is SPM.xBF.T0 = 8 (middle slice), but if you choose the first slice
			  %	       as reference you have to change SPM.xBF.T0 = 1. 
			  %	       The parameter can be find in the parameter-block for defining the design matrix.
 
dir_base        = '/Users/pliu/Documents/LayerPRF/LayerPRF'; % directory which contains all subject-specific directory
dir_analysis    = 'Results';   % directory where SPM.mat, beta images, con images etc. are stored for each subject
dir_functional  = 'Functional/sess1'; % directory of functional data (NIfTI files)
dir_anatomical  = 'Anatomical'; % directory including the T1 3D high resoluation anatomical data set
name_scans      = 'data';       % name of the NIfTI functional data-file
name_fmap	    = 'fmap';       % name of the NIfTI fieldmap-file
name_ana        = 'T1';        % name of the (high resolution) anatomical NIfTI file 
sess_prfx       = 'sess';       % prefix of session directories (only needed for multi-session analyses)

delete_files	= 0;  % delete intermediate steps during pre-processing? 1=yes, 0=no
data_scale      = 1; % scale data before realign and unwarp: 
                      % data scaling factor, ENSURE THAN MAXUMUM VALUE WILL BE < 32000
unwarp          = 1;  % unwarp during realignment? 1=yes, 0=no
fieldmap 	    = 0;  % fieldmap scan available? 1=yes, 0=no
slice_timing    = 2;  % slice timing? 2=before realignment, 1=after realignment, 0=no correction  (because interleaved acquisition)
slice_order     = 2;  % slice order? 2=custom, 1=ascending, 0=descending
slice_ord       = [2:2:30 1:2:29];	% Syntax: [first_slice:increment/decrement:last_slice]
							% Examples:
							% descending acquisition [num_slices:-1:1]
							% ascending acquisition [1:1:num_slices]
							% interleaved acquisition (ascending - even slices first): [2:2:num_slices 1:2:num_slices-1]
norm_type       = 1;  % normalization: 1=T1 segmentation, 0=EPI template
norm_preorient  = 1;  % pre-orient all images with the MNI space by co-registration with T1 template
dst_reso_func   = 1;% resolution for functional data after normalization in mm
dst_reso_ana    = 1;% resolution for anatomical data after normalization in mm
smooth_fwhm     = 0;  % smoothing after normalization in mm (at minimum twice the voxel size)
rp              = 0;  % include realignment parameters in design matrix? 1=yes, 0=no
cutoff_highpass = 128;% cutoff for high pass filter in s [Inf = no filtering]
start_analysis  = 3;  % start of analysis
    % 1 = realignment & unwarp and slicetime
    % 2 = normalization & smoothing
    % 3 = design matrix definition & parameter estimation
    % 4 = define and compute contrasts
    
% names of subjects (list all names of subject directories)
%name_subj = {'fip66_4571/blocked_design1','fip66_4571/blocked_design2','mhc403_4574/blocked_design1','mhc403_4574/blocked_design2'};
name_subj = {'rrl733'};

% definition of conditions
cond_names   = {'D2','D3'};
cond_param   = [ 0 0 ];    % number of parameters per condition
% NB: prepare subject+session specific onset and duration files 'VISUAL.ons', 'VISUAL.dur', 'AUDIO.ons', 'AUDIO.dur' 
% and put this files in the *same* directory where you have the raw data 'data.nii'. 
% If your condition is an event just put a single 0 into the duration file.
% If there are parameters for conditions, specify parameter files (e.g. 'AUDIO.par').
% ALL vectors must be written as columns.

% definition of contrasts
c_name = {'D2' ; 'D3'};
c_con  = [ 1 -1 ; -1 1 ];
c_type = {'T','T'};


%===========================================================================
% end of user specified settings
%===========================================================================





%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
%
%       Do not edit below.
%
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------


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
    dir_ana = [dir_base fs name_subj{k} fs dir_anatomical];
    
    
    % realign & unwarp
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
       
   	    % load NIfTI files
   	    %-----------------------------
   	    cd (dir_scans{sess})
   	    P = fullfile(dir_scans{sess},[prefix name_scans '.nii']);
	    
    	
    	    % Realignment (coregister only)
    	    %------------------------------------
    	    FlagsC  = struct('quality',1,'fwhm',5,'sep',4,'interp',2,'wrap',[0 0 0],'rtm',0,...
	                     'PW','','graphics',1,'lkp',1:6);
    	    spm_realign(P,FlagsC);
    	
    	    % Reslice (unwarp and reslice)
    	    %------------------------------------
    	    if ~unwarp
	
    	    	FlagsR  = struct('interp',1,'mask',1,'mean',1,'which',2,'wrap',[0 0 0]','prefix','r');
    	    	spm_reslice(P,FlagsR);		 
    	
    	    else 

		if fieldmap == 1	% fieldmap scan available?
    			fmap = fullfile(dir_scans{sess},[name_fmap '.nii']);			
		else
                        fmap = [];
        	end

    		% Estimate unwarping parameters
		uw_est_flags = struct('order',[12 12],'sfP',fmap,'regorder',1,'jm',0,'fot',[4 5],...
		                      'sot',[],'fwhm',4,'rem',1,'exp_round','Average','noi',5,'hold',[1 1 1 0 1 0]);
    		ds     = spm_uw_estimate(P,uw_est_flags);
    		pefile = fullfile(dir_scans{sess},[name_scans '_uw.mat']);
    		save(pefile,'ds');
    		    
    		% Write unwarped images
    		uw_write_flags = struct('mask',1,'mean',1,'interp',4,'wrap',[0 1 0],'which',1,'udc',1,'prefix','u');		
    		spm_uw_apply(ds,uw_write_flags);    
    		
    	    end % unwarp
    
    	    % Slice-timing correction 
    	    %------------------------
    	    if slice_timing == 1
    		   		
    		% load NIfTI files
    		%-----------------------------
		if unwarp
    		   prefix = 'u';
		else
		   prefix = '';
		end 
		cd (dir_scans{sess})
		P = fullfile(dir_scans{sess},[prefix name_scans '.nii']);     	
     	 
    		% Slicetime correction
    		%------------------------------------
    		spm_slice_timing(P, slice_ord, ref_slice, [slice_time slice_gap], 'a');

    	    end % slice timing

    	    % delete intermediate steps (I)
    	    %------------------------------
    	    if delete_files == 1
		if (slice_timing > 0 && unwarp == 1)
    		    if slice_timing == 2
	    		del_prefix = 'a';
		    elseif slice_timing == 1   
			if unwarp
	    		    del_prefix = 'u';
	    		else
	    		    del_prefix = 'r';
	    		end
		    end
		del = fullfile(dir_scans{sess},[del_prefix name_scans '.*']);
	        delete(del);
		end
	    end
	
	end % loop over sessions
    
    end % if-block for realign & unwarp section
    

    % Normalization & smoothing
    %=======================================
    
    if start_analysis <= 2
    
	for sess = 1:n_sess,
            	    	
    	    % Load NIfTI files
    	    %-----------------------------
	    if unwarp
    	       prefix = 'u';
	       meanprefix = 'u';
	    else
	       prefix = '';
	       meanprefix = '';
	    end
	    if slice_timing==1
	       prefix = ['a' prefix];
	    end
	    if slice_timing==2
	       prefix = [prefix 'a'];
	       meanprefix = [meanprefix 'a'];
	    end	   
	    cd (dir_scans{sess})
	    mean_img = fullfile(dir_scans{sess},['mean' meanprefix name_scans '.nii']);
	 
    	    % Estimate normalisation parameters
    	    %----------------------------------
	    if norm_type==0

         	    sn_estimate_flags = struct('smosrc',8,'smoref',0,'regtype','mni','cutoff',30,...
		                               'nits',16,'reg',0.1,'graphics',1);
	    	    matname = [spm_str_manip(mean_img,'sd') '_sn.mat'];
    		    VG      = fullfile(spm('Dir'),'templates',['EPI' '.nii']);
		    sn      = spm_normalise(VG,mean_img,matname,'','',sn_estimate_flags);

	    else
	    
		    % Coreg
		    %-------------------------
		    copyfile(fullfile(dir_ana,[name_ana '_pl.nii']),dir_scans{sess});	    
		    copyfile(fullfile(dir_ana,[name_ana '_np.nii']),dir_scans{sess});	    
		    copyfile(fullfile(dir_ana,[name_ana '_7T.nii']),dir_scans{sess});	    

		    struc_7T = fullfile(dir_scans{sess},[name_ana '_7T.nii']);
		    struc_pl = fullfile(dir_scans{sess},[name_ana '_pl.nii']);
		    struc    = fullfile(dir_scans{sess},[name_ana '_np.nii']);


                    coreg_flags = struct('sep',[4 2],'params',[0 0 0  0 0 0],... 
		    			 'cost_fun','ncc','fwhm',[7 7],'tol',...
			   		 [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001],...
			   		 'graphics',1);	

                    % coreg peeld 3T image onto peeled 7T image
		    %---------------------------------------------
		    x  = spm_coreg(struc_7T, struc_pl, coreg_flags);		 
		    M  = inv(spm_matrix(x));
		    
		    MM = spm_get_space(deblank(struc_pl));			 
		    spm_get_space(deblank(struc_pl), M*MM);		 

		    MM = spm_get_space(deblank(struc));			 
		    spm_get_space(deblank(struc), M*MM);		 
					 
					 					 	    		    
		    %coreg peeled 3T image onto mean image
                    %-------------------------------------
		    x  = spm_coreg(mean_img, struc_pl);
		    M  = inv(spm_matrix(x));
		    
		    MM = spm_get_space(deblank(struc_pl));
		    spm_get_space(deblank(struc_pl), M*MM);
		    	
		    MM = spm_get_space(deblank(struc));
		    spm_get_space(deblank(struc), M*MM);
		    	
		    MM = spm_get_space(deblank(struc_7T));
		    spm_get_space(deblank(struc_7T), M*MM);
			
			
		    % pre-orient all images with the MNI space by co-registration with T1 template
		    %-----------------------------------------------------------------------------	
		    if norm_preorient==1	
                         coreg_flags = struct('sep',[4 2],'params',[0 0 0  0 0 0],... 
			                      'cost_fun','ncc','fwhm',[7 7],'tol',...
					      [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001],...
					      'graphics',1);
			 T1_img = fullfile(spm('dir'),'templates','T1.nii');		      
			 x  = spm_coreg(T1_img, struc, coreg_flags);
    	                 M  = inv(spm_matrix(x));

			 % struct
			 %-------
    	                 MM = spm_get_space(deblank(struc));
    	                 spm_get_space(deblank(struc), M*MM);

    	                 MM = spm_get_space(deblank(struc_pl));
    	                 spm_get_space(deblank(struc_pl), M*MM);

    	                 MM = spm_get_space(deblank(struc_7T));
    	                 spm_get_space(deblank(struc_7T), M*MM);
			 
                         % mean img
			 %---------
    	                 MM = spm_get_space(deblank(mean_img));
    	                 spm_get_space(deblank(mean_img), M*MM);			 
			 
                         % timeseries
			 %-----------
			 for i=1:n_scans
			    P = fullfile(dir_scans{sess},[prefix name_scans '.nii,' mat2str(i)]);
			    MM = spm_get_space(deblank(P));
    	                    spm_get_space(deblank(P), M*MM);
			 end
		    end	
				    
		    % 1st pass
                    % -------------------------------
                    estopts.regtype=''; 		    % turn off affine
                    out = spm_preproc(struc,estopts);
                    [sn,isn]  = spm_prep2sn(out);

                    % only write out attenuation corrected image
                    writeopts.biascor = 1;
                    writeopts.GM  = [0 0 0];
                    writeopts.WM  = [0 0 0];
                    writeopts.CSF = [0 0 0];
                    writeopts.cleanup = [0];
                    spm_preproc_write(sn,writeopts);
                    
                    [pth,nam,ext,num] = spm_fileparts(struc);
                    struc = fullfile(pth,['m' nam ext]);	    


                    % 2nd pass
                    % -------------------------------
                    estopts.regtype='mni';		    % turn on affine
                    out = spm_preproc(struc,estopts);
                    [sn,isn]  = spm_prep2sn(out);

                    % assume GM(2) means unmod
                    writeopts.biascor = 1;
                    writeopts.GM  = [0 1 1];
                    writeopts.WM  = [0 0 0];
                    writeopts.CSF = [0 0 0];
                    writeopts.cleanup = [0];
                    spm_preproc_write(sn,writeopts);

                    save(sprintf('%s_seg_sn.mat',spm_str_manip(struc,'sd')),'sn')
                    save(sprintf('%s_seg_inv_sn.mat',spm_str_manip(struc,'sd')),'isn')
		    
		    [pth,nam,ext,num] = spm_fileparts(struc);
                    struc = fullfile(pth,['m' nam ext]);
                     
	    end	    

	    
    	    % Write normalised & smoothed
    	    %-------------------------------------------------
	    sn_write_flags = struct('interp',1,'vox',[dst_reso_func dst_reso_func dst_reso_func],...
	                            'bb',[-70 -108 15 ; 70 30 75],'wrap',[0 0 0],'preserve',0,'prefix','w');


 % [-78 -112 20 ; 78 76 85]
 % [-78 -112 -50 ; 78 76 85]
% [-70 -108 15 ; 70 68 75]
            clear P;
            P = fullfile(dir_scans{sess},[prefix name_scans '.nii']);
    	    msk = spm_write_sn(P,sn,sn_write_flags,'mask');
    	    spm_write_sn(mean_img,sn,sn_write_flags,msk);
    	    spm_write_sn(P,sn,sn_write_flags,msk);
	    if norm_type==1	    
	       sn_write_flags_ana = struct('interp',1,'vox',[dst_reso_ana dst_reso_ana dst_reso_ana],...
	                                   'bb',[-78 -112 -50 ; 78 76 85],'wrap',[0 0 0],'preserve',0,'prefix','w');		  	    	          
               spm_write_sn(struc,sn,sn_write_flags_ana);
	       spm_write_sn(struc_pl,sn,sn_write_flags_ana);
	       spm_write_sn(struc_7T,sn,sn_write_flags_ana);
	    end

    	    % delete intermediate steps (II)
    	    %-------------------------------
    	    if delete_files == 1
    		if unwarp
    	    	   del_prefix = 'u';
		else
	    	   del_prefix = '';
		end
		if slice_timing==1
	           del_prefix = ['a' del_prefix];
	        end
		if slice_timing==2
	           del_prefix = [del_prefix 'a'];
		end	   
		del = fullfile(dir_scans{sess},[del_prefix name_scans '.*']);
		delete(del);
	    end

    	    % smooth images
    	    %------------------------------------------------
	    cd (dir_scans{sess})
    	    prefix = ['w' prefix];
	    
	    fname = fullfile(dir_scans{sess},['s' prefix name_scans '.nii']);
	    P = fullfile(dir_scans{sess},[prefix name_scans '.nii']);
            spm_smooth(P,fname,smooth_fwhm); 
	    
    	    % delete intermediate steps (III)
    	    %--------------------------------
    	    if delete_files == 1
	        del_prefix = ['w' del_prefix];
		del = fullfile(dir_scans{sess},[del_prefix name_scans '.*']);
		delete(del);
	    end

	end % loop over sessions	    
    		
    end    % if-block for normalisation & smoothing
   
     
    % define design matrix
    %===========================================================================
    if start_analysis <= 3

        cd ([dir_base fs name_subj{k} fs dir_analysis]);
        
        % number of scans, sessions, TR, basis functions and timing parameters
        %---------------------------------------------------------------------------
	SPM.nscan          = ones(1,n_sess)*n_scans;
        SPM.xY.RT          = TR;
        SPM.xBF.name       = 'hrf';		% basis functions and timing parameters
						% OPTIONS:'hrf'
						%         'hrf (with time derivative)'
						%         'hrf (with time and dispersion derivatives)'
						%         'Fourier set'
						%         'Fourier set (Hanning)'
						%         'Gamma functions'
						%         'Finite Impulse Response'
        SPM.xBF.order      = 1;			% order of basis-set              
        SPM.xBF.length     = 32;                % length in seconds
        SPM.xBF.T          = 16;           	% number of time bins per scan
	SPM.xBF.dt         = TR/SPM.xBF.T;
        SPM.xBF.T0         = 8;			% middle slice/timebin
        SPM.xBF.UNITS      = 'secs';            % OPTIONS: 'scans'|'secs' for onsets
        SPM.xBF.Volterra   = 1;                 % OPTIONS: 1|2 = order of convolution
        SPM.xBF            = spm_get_bf(SPM.xBF);

        % Trial specification: onsets, duration (in scans) and 
        % parameters for modulation
        %---------------------------------------------------------------------------
        for sess=1:n_sess
        	for c=1:length(cond_names)
         		SPM.Sess(sess).U(c).name{1}   = cond_names{c};
         		SPM.Sess(sess).U(c).ons       = (load (fullfile(dir_scans{sess},[cond_names{c} '.ons'])));
         		SPM.Sess(sess).U(c).dur       = (load (fullfile(dir_scans{sess},[cond_names{c} '.dur'])));
			SPM.Sess(sess).U(c).P(1).name = 'none';
			SPM.Sess(sess).U(c).P(1).h    = 0;
			if cond_param(c)>0
			   for param = 1:cond_param(c)
			      AllParam = load (fullfile(dir_scans{sess},[cond_names{c} '.par']));
			      SPM.Sess(sess).U(c).P(param).P    = AllParam(:,param);
			      SPM.Sess(sess).U(c).P(param).name = ['param' num2str(param)]; 
			      SPM.Sess(sess).U(c).P(param).h    = 1;
			      SPM.Sess(sess).U(c).name{param+1} = ['param' num2str(param)];
			   end
			end   
        	end
        end
        
        % User-specified regressors
        %--------------------------
        if rp
            % realignment parameters
            rnames = {'translX','translY','translZ','rotX','rotY','rotZ'};
            for sess = 1:n_sess
	        if slice_timing==2
		   fn  = fullfile(dir_scans{sess},['rp_a' name_scans '.txt']);
		else
		   fn  = fullfile(dir_scans{sess},['rp_' name_scans '.txt']);
		end      
                [r1,r2,r3,r4,r5,r6] = textread(fn,'%f%f%f%f%f%f');
                cd ([dir_base fs name_subj{k} fs dir_analysis]);
                SPM.Sess(sess).C.C = [r1 r2 r3 r4 r5 r6];
                SPM.Sess(sess).C.name = rnames;
            end
        else
            % no user specified regressors 
            for sess=1:n_sess
                SPM.Sess(sess).C.C = [];
                SPM.Sess(sess).C.name = {};
            end
        end

        % Define design matrix, global normalization, high-pass and intrinsic correlation
        %---------------------------------------------------------------------------
        [SPM] = spm_fMRI_design(SPM);
        save SPM_DesignMatrix SPM
        SPM.xGX.iGXcalc    = 'None';		% global normalization: OPTIONS:'Scaling'|'None'
	SPM.xVi.form	   = 'AR(1)';		% intrinsic autocorrelations: OPTIONS: 'none'|'AR(1)'
        for sess=1:n_sess
            SPM.xX.K(sess).HParam = cutoff_highpass;
        end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% prefix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%----------
	%if unwarp
	 %   prefix = 'u';
	%else
	 %   prefix = 'r';
	%end		    
	%if slice_timing==1
	 %   prefix = ['a' prefix]; 
	%end
	%if slice_timing==2
	 %   prefix = [prefix 'a']; 
	%end	       
	prefix = 'a';  

        % specify data by matrix of filenames and configure design matrix
        %------------------------------------------------------------------------
	P	    = [];
        for sess    = 1:n_sess,
	    cd (dir_scans{sess})
	    tmpP = fullfile(dir_scans{sess},[prefix name_scans '.nii']);
	    P = [P ; tmpP];
        end		  
        cd ([dir_base fs name_subj{k} fs dir_analysis]);
	SPM.xY.P    = P;

        % Configure design matrix and estimate parameters
        %------------------------------------------------
        SPM = spm_fmri_spm_ui(SPM);
	SPM = spm_spm(SPM);

    end % if-block for design matrix definition & parameter estimation


    
    % Add contrasts
    %===========================================================================
    if start_analysis <= 4
                
        % switch to analysis directory and load SPM
        %------------------------------------------
        clear SPM;
        cd ([dir_base fs name_subj{k} fs dir_analysis]);
        load SPM;	     

        % Create 1st contrast for 'effects of interest'
        %-----------------------------------------------
        try
        	iX0     = [SPM.xX.iB SPM.xX.iG];
        catch
        	iX0     = [];
        end
        SPM.xCon        = spm_FcUtil('Set','effects of interest','F','iX0',iX0,SPM.xX.xKXs);

        % Add further contrasts
        %------------------------------------------------------------------
        n_cntr      = length(SPM.xCon); 			% number of contrasts defined so far
        eb          = eye(SPM.xBF.order);			% basis functions     
        for c = 1:length(c_name)
	    c_wgt{c}           = [kron(c_con(c,:),eb(1,:))];
	    c_complete	       = [];
	    for sess = 1:n_sess
	        if rp
	           c_complete  = [c_wgt{c} 0 0 0 0 0 0 c_complete];
		else
		   c_complete  = [c_wgt{c} c_complete];   
		end   
	    end
            cw                 = [c_complete zeros(size(c_wgt{c},1), n_sess)]'	% pad with zeros for constants
            SPM.xCon(c+n_cntr) = spm_FcUtil('Set',c_name{c},c_type{c},'c',cw,SPM.xX.xKXs);
        end
        
        % Evaluate contrasts
        %---------------------------------------------------------------------------
        spm_contrasts(SPM);

    end     % if-block for contrast definition & evaluation
    
    
    
    % Switch to next subject
    %=======================================
    k   = k + 1;
    clear SPM P tmpP;
    
    
end     % end of main loop

cd (dir_base);
