% GET NAMES OF THE SUBFOLDERS IN A FOLDER

% COPYRIGHT: 
%
%   MIT License
%   © 2023 Dario Chiantello <pipesoundlibrary@gmail.com>
%

function folderNames = getFolderNames(path)

    if exist(path, 'dir')
      
        folderFiles                     = dir(path);
        dirFlags                        = [folderFiles.isdir];
        folders                         = folderFiles(dirFlags);

        folderNames                     = cell(length(folders)-2,1);
        [folderNames{:}]                = deal(folders(3:end).name);
        
    else
        
        folderNames = [];
        
    end

