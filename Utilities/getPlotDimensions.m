% PLOT DIMENSIONS

% COPYRIGHT: 
%
%   MIT License
%   © 2023 Dario Chiantello <pipesoundlibrary@gmail.com>
%

function size = getPlotDimensions(propertyString)


    % DEFAULT 
    
    %default
    if strcmp(propertyString, 'default_figurePosition')

        size       = [0, 0, 24, 12];   

    elseif  strcmp(propertyString, 'default_paperSize')
        
        size       = [24, 13];   
        
        

        
    % CUSTOM
        
    %custom_full_1
    elseif strcmp(propertyString, '165x52_figurePosition')
    
        size       = [0, 0, 16.5, 5.2]; 

    elseif strcmp(propertyString, '165x52_paperSize')

        size       = [16.5, 5.2]; 
        
        
    %custom_threeQuarters_1
    elseif strcmp(propertyString, '125x67_figurePosition')

        size       = [0, 0, 12.5, 6.7]; 

    elseif strcmp(propertyString, '125x67_paperSize')

        size       = [12.5, 6.7]; 
    
        
    %custom_half_1
    elseif strcmp(propertyString, '65x45_figurePosition')

        size       = [0, 0, 6.5, 4.5];

    elseif strcmp(propertyString, '65x45_paperSize')

        size       = [6.5 4.5];  
        
        
         
              
    %IEEE
    
    %ieee_full_1
    elseif strcmp(propertyString, '181x565_figurePosition')

        size       = [0, 0, 18.1, 5.6];

    elseif strcmp(propertyString, '181x56_paperSize')

        size       = [18.1, 5.6];  
    
    %ieee_threeQuarters_1
    elseif strcmp(propertyString, '160x56_figurePosition')

        size       = [0, 0, 14, 5.6];

    elseif strcmp(propertyString, '160x56_paperSize')

        size       = [14, 5.6];  
            
    %ieee_half_1
    elseif strcmp(propertyString, '88x50_figurePosition')

        size       = [0, 0, 8.8, 5];

    elseif strcmp(propertyString, '88x50_paperSize')

        size       = [8.8, 5];  
         

    end


