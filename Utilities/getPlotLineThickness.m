% PLOT LINE THICKNESS

% COPYRIGHT: 
%
%   MIT License
%   © 2023 Dario Chiantello <pipesoundlibrary@gmail.com>
%

function thickness = getPlotLineThickness(propertyString)

    if     strcmp(propertyString, 'thinLine') == true

        thickness       = 0.5;

    elseif strcmp(propertyString, 'thickLine') == true

        thickness       = 1.2;

    else
        
        thickness       = 0.5;
        
    end
    
end


