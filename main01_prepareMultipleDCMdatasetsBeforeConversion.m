clear all;
clc;

rootDir = uigetdir(pwd, 'Select DCM directory...');
if rootDir == 0
    return 
end

outputDirForFakeRefImages = uigetdir(pwd, 'Select output directory for fake ref. images...');
if outputDirForFakeRefImages == 0
    return 
end

outputDirForMoCoFiles = uigetdir(pwd, 'Select output directory for MoCoFiles...');
if outputDirForMoCoFiles == 0
    return 
end


disp('LIST OF DIRECTORIES:');
dirList = getListOfDirectoriesWithinRootDir(rootDir);
disp(' ');

for dirIdx = 1:length(dirList)
    disp(['process folder -> ' dirList(dirIdx).path]);
    prepareOneSingleDCMdatasetBeforeConversion(dirList(dirIdx).path, outputDirForFakeRefImages, outputDirForMoCoFiles);
    disp(' ');
end

