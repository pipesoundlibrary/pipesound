% SET SCALE FACTOR AND UNIT IN THE AUDIO FILE METADATA

% COPYRIGHT: 
%
%   MIT License
%   © 2023 Dario Chiantello <pipesoundlibrary@gmail.com>
%

function commentString = setUnitAndScale(scalefactor, unit)


commentString = sprintf('scaleFactor%s=%.6f', unit, scalefactor);