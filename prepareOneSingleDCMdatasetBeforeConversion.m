function prepareOneSingleDCMdatasetBeforeConversion(dicomDir, fakeRefImgDir, MoCoFilesDir)

filelist = dir(dicomDir);

% print patient ID that belongs to this folder (just for one file)
for fileIdx = 1:length(filelist)
    if (filelist(fileIdx).name(1) ~= '.') && (filelist(fileIdx).isdir == 0) && (~strcmpi(filelist(fileIdx).name, 'DICOMDIR')) && (strcmp(filelist(fileIdx).name(1:3), 'MR.'))
        
        % read out all properties of the current file
        dcmInfo = dicominfo([dicomDir filesep filelist(fileIdx).name]);

        % disp patient ID
        if isfield(dcmInfo, 'PatientID')
            disp(['   ***dcmInfo.PatientID: ' dcmInfo.PatientID]);
            break
        end
    end
end



for i = 1:1000
    series(i).mocoParam = [];
end

for fileIdx = 1:length(filelist)
    
    
%     if (filelist(fileIdx).name(1) ~= '.') && (filelist(fileIdx).isdir == 0) && ~strcmpi(filelist(fileIdx).name, 'DICOMDIR')
    if (filelist(fileIdx).name(1) ~= '.') && (filelist(fileIdx).isdir == 0) && ~strcmpi(filelist(fileIdx).name, 'DICOMDIR') && (strcmpi(filelist(fileIdx).name(1:3), 'MR.') || strcmpi(filelist(fileIdx).name(1:3), 'SRe'))
       
        % read out all properties of the current file
        dcmInfo = dicominfo([dicomDir filesep filelist(fileIdx).name]);

        % Look out for Reference Volumes which have to be excluded
        if isfield(dcmInfo, 'ImageComments') && ((strcmpi(dcmInfo.ImageComments, 'Fake MoCo Reference Volume') || (strcmpi(dcmInfo.ImageComments, 'Reference volume for motion correction. + DiCo Applied'))))
            
            % move fake ref image to the designated directory
            movefile([dicomDir filesep filelist(fileIdx).name], fakeRefImgDir)
            disp([dcmInfo.ImageComments ' -> ' filelist(fileIdx).name]);        
            
        % Look out for MoCo-Images, of which we want to extract the
        % movement parameters
        elseif strcmpi(dcmInfo.SeriesDescription, 'MoCoSeries_DiCo')
            curFile_seriesNr = dcmInfo.SeriesNumber;
            curFile_patientID = dcmInfo.PatientID;
            curFile_acqNr = dcmInfo.AcquisitionNumber;

            % Save MoCo info in a variable that can be used by SPM as additional regressors
            C = strsplit(dcmInfo.ImageComments, {',', ' '});
            
            % the moco parameters in now stored in the order as the dcm-files
            % are processed (the order is determined by the filename).
            % but sometimes the order of the files (resp. filenames) is not identical to the
            % acquisition order of the volumes. So we want to sort the parameters later by
            % making use of the acqisition number that is stored in the
            % DCM-header info. For that reason, we store the acquisition
            % number too (in column 1).
            series(curFile_seriesNr).mocoParam(end+1, :) = [curFile_acqNr, str2double(C{2}), str2double(C{3}), str2double(C{4}), str2double(C{5}), str2double(C{6}), str2double(C{7})];

        end
        
    end
end

% Save MoCo info in a mat-file that can be used by SPM as additional regressors
for i = 1:length(series)
    if ~isempty(series(i).mocoParam)
        
        R = series(i).mocoParam;
        R = sortrows(R,1); % sort the parameters in the order of the acquisition number (column 1)
        R = R(:,2:end); % delete column 1 (acquisition number), it is not used anymore.
        
        save([MoCoFilesDir filesep curFile_patientID '_MoCoParam_7T_Series' num2str(i) '.mat'], 'R');
    end
end


