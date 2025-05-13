%%% this script requires SPM toolbox
%%% user specified path to the SPM8 directory
%%% this script converts DICOM files to NIFTI


clear all;
clc;

rootDir = uigetdir(pwd, 'Select input direct');
if rootDir == 0
    return 
end

% select directory where DCM files a
outputDir = uigetdir(pwd, 'Select output directory (for *.nii files)...');
if outputDir == 0
    return 
end


disp('LIST OF DIRECTORIES:');
dirList = getListOfDirectoriesWithinRootDir(rootDir);
disp(' ');

for dirIdx = 1:length(dirList)
    disp(['process folder -> ' dirList(dirIdx).path]);
    
    listOfRawScans = dir([dirList(dirIdx).path filesep 'MR.*']);
    
    
    if ~isempty(listOfRawScans)
        
        disp(['   convert ' num2str(length(listOfRawScans)) ' DCM files to *.nii']);
        
        % BRING THE LIST OF DCM FILES IN THE SPM FILELIST-FORMAT
        for i = 1:length(listOfRawScans)
            spmFilelist{i,1} = [dirList(dirIdx).path filesep listOfRawScans(i).name];
        end

        % DEFINE SPM JOB
        matlabbatch{1}.spm.util.dicom.data = spmFilelist;
      
        matlabbatch{1}.spm.util.dicom.root = 'patid';
        matlabbatch{1}.spm.util.dicom.outdir = {[outputDir filesep]};
        matlabbatch{1}.spm.util.dicom.convopts.format = 'nii';
        matlabbatch{1}.spm.util.dicom.convopts.icedims = 0;

        % RUN SPM JOB
        spm_jobman('initcfg');
        spm_jobman('run', matlabbatch);

        disp(' ');
    end
end
