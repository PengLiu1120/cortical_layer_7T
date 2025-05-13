%%% this script requires the imaging processing toolbox from matlab
%%% this script sorts the directories according to the variables
%%% "ProtocolName" and "SeriesDescription" of the DICOM header

%%% note: two field mapping series here in one folder

%%% written by Esther Kuehn 22 January 2019
%%% contact: esther.kuehn@dzne.de

% select directory where DCM files are
% this will also the the directory where the new folders are generated
dicomDir = uigetdir(pwd, 'Select DCM directory...');
if dicomDir == 0
    return 
end

FilesDir = dicomDir;

% first loop: create folders and copy data into folders
filelist = dir(dicomDir);

 for fileIdx = 1:length(filelist)
     if (filelist(fileIdx).name(1) ~= '.') && (filelist(fileIdx).isdir == 0) && (~strcmpi(filelist(fileIdx).name, 'DICOMDIR')) && (strcmp(filelist(fileIdx).name(1:3), 'MR.'))
        
         % read out all properties of the current file
        dcmInfo = dicominfo([dicomDir filesep filelist(fileIdx).name]);
        
        % Create folders for each Series Description if it does not exist
        % already
        mkdir(dicomDir,[dcmInfo.ProtocolName '_' dcmInfo.SeriesDescription]); 
        
        % Move file to folder
        folder = [dicomDir '/' dcmInfo.ProtocolName '_' dcmInfo.SeriesDescription '/'];
        movefile(dcmInfo.Filename,folder);
         
     end 
 end


 % create folder called PhoenixZIPReport and move images there that start
 % with "SRe"    
  mkdir(dicomDir,['PhoenixZIPReport']);
  filelist2 = dir(dicomDir);
  
  for fileIdx = 1:length(filelist2)
     if (filelist2(fileIdx).name(1) ~= '.')  && (filelist2(fileIdx).isdir == 0) && (strcmp(filelist2(fileIdx).name(1:3), 'SRe'))
     
         % read out all properties of the current file
        dcmInfo3 = dicominfo([dicomDir filesep filelist2(fileIdx).name]);
           
        % Move file to folder
        folder = [dicomDir '/PhoenixZIPReport'];
        movefile(dcmInfo3.Filename,folder);
         
     end 
 end
