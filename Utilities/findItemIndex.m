% FIND THE ITEM INDEX IN A LIST BY SPECIFYING ITS NAME

% COPYRIGHT: 
%
%   MIT License
%   © 2023 Dario Chiantello <pipesoundlibrary@gmail.com>
%

function itemIndex = findItemIndex(list, name)

    itemIndex = find(cellfun(@(x) strcmp(x,name), list), 1);

end