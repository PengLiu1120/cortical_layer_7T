function dirList = addNewSubdirsToDirList(currentDirToCheck, dirList)

dirInfo = dir(currentDirToCheck);

for i = 1:length(dirInfo)
    if(dirInfo(i).isdir == 1) && (dirInfo(i).name(1)~='.')
        dirIsAlreadyInList = false;
        for dirListElem = 1:length(dirList)
            if strcmpi([currentDirToCheck filesep dirInfo(i).name], dirList(dirListElem).path)
                dirIsAlreadyInList = true;
            end
        end
        if ~dirIsAlreadyInList
            dirList(end+1).path = [currentDirToCheck filesep dirInfo(i).name];
%             disp(dirList(end).path);
        end
    end
end
