% GET COLOUR BY INDEX

% COPYRIGHT: 
%
%   MIT License
%   © 2023 Dario Chiantello <pipesoundlibrary@gmail.com>
%

function colour = getColourByIndex(index)

if index == 0

    colour = getColour('blue'); 
    
elseif index == 1
        
    colour = getColour('green'); 
        
elseif index == 2
        
    colour = getColour('lightGreen'); 
    
elseif index == 3
        
    colour = getColour('red'); 
    
elseif index == 4
        
    colour = getColour('azure'); 
    
elseif index == 5
        
    colour = getColour('yellow'); 
    
elseif index == 6
        
    colour = getColour('orange'); 
    
elseif index == 7
        
    colour = getColour('violet'); 
    
elseif index == 8
        
    colour = getColour('granate');  
    
elseif index == 9
        
    colour = getColour('gray'); 
    
elseif index == 10
        
    colour = getColour('pink'); 
    
elseif index == 11
        
    colour = getColour('peach'); 
    
elseif index == 12
        
    colour = getColour('darkViolet'); 

else
    
    colour = getColour('black'); 
    
end 