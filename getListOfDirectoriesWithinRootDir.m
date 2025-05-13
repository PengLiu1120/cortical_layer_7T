function dirList = getListOfDirectoriesWithinRootDir(rootDir)
% function creates list of all subdirectories, subsubdirectories and so on
% within a given root directory (rootDir)
% all directories are listed in the output variable dirList, which you can
% access via dirList(i).path
% The root directory is also in the output variable dirList (field 1).

dirList(1).path = rootDir;

newDirFound = true;
newDirList = dirList;

while newDirFound
    
    for dirElem = 1:length(dirList)
        currentDirToCheck = dirList(dirElem).path;
        newDirList = addNewSubdirsToDirList(currentDirToCheck, newDirList);
    end
    
    if length(dirList) == length(newDirList)
        newDirFound = false;
    end
    
    dirList = newDirList;
    
end

for i = 1:length(dirList)
    disp(['   -> ' dirList(i).path]);
end
