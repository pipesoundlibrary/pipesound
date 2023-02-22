% FILE WRITE AND PROMPT

% COPYRIGHT: 
%
%   MIT License
%   © 2023 Dario Chiantello <pipesoundlibrary@gmail.com>
%

function printFileAndPrompt(fileId, outputString, colourString)

    if ~isempty(fileId)
    
        fprintf(fileId, outputString);
    
    end
    
    if ~exist('colourString', 'var')
        
        fprintf(outputString);
   
    else
    
        cprintf(['*', mat2str(getColour(colourString))], outputString);
        
    end
end
