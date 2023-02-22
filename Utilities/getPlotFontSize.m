% FONT SIZE

% COPYRIGHT: 
%
%   MIT License
%   © 2023 Dario Chiantello <pipesoundlibrary@gmail.com>
%

function fontSize = getPlotFontSize(propertyString)

    if     strcmp(propertyString, 'defaultLegendFontsize')

        fontSize       = 10;
    
    elseif strcmp(propertyString, 'defaultLabelFontsize')

        fontSize       = 10;
             
    elseif strcmp(propertyString, 'defaultTickFontsize')
       
        fontSize       = 10;
        
        
    elseif strcmp(propertyString, 'ieeeLegendFontSize')    
        
        fontSize       = 8;
        
    elseif strcmp(propertyString, 'ieeeLabelFontSize')

     	fontSize       = 8;
        
    elseif strcmp(propertyString, 'ieeeTickFontSize') 
        
        fontSize       = 8;
        
        
    else
        
        fontSize       = 13;
        
    end
    
end