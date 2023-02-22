% GET SCALE FACTOR AND UNIT FROM THE AUDIO FILE METADATA

% COPYRIGHT: 
%
%   MIT License
%   © 2023 Dario Chiantello <pipesoundlibrary@gmail.com>
%


function [scaleFactor, unit] = getUnitAndScale(metadata)

    comment = metadata.Comment;

    if isempty(comment)
        
        scaleFactor = 1;
        warning('No scale factor and no physical unit found')
        
    else

        try

            i = strfind(comment, '=');
            scaleFactor = metadata.Comment(i+1:end);
            scaleFactor = str2double(scaleFactor);

        catch

            scaleFactor = 1;
            warning('No valid scale factor found')

        end

        
        try
    
            b1 = strfind(comment, '[');
            b2 = strfind(comment, ']');

            unit = comment(b1:b2);

            if isempty(unit)

                unit = '';
                warning('No physical unit found')

            end

        catch

            unit = '';
            warning('No physical unit found')

        end
        
    end
    
    
    
   