% GET COLOURS

% COPYRIGHT: 
%
%   MIT License
%   © 2023 Dario Chiantello <pipesoundlibrary@gmail.com>
%

function colour = getColour(colourString)


if strcmp(colourString, 'black')

    colour       = [0,      0,      0];
    
elseif strcmp(colourString, 'green')
        
    colour       = [0,      0.4,    0.3];
        
elseif strcmp(colourString, 'lightGreen')
        
    colour       = [0.4660, 0.6740, 0.1880]; 
    
elseif strcmp(colourString, 'red')
        
    colour       = [0.8500, 0.3250, 0.0980];
    
elseif strcmp(colourString, 'blue')
        
    colour       = [0,      0.4470, 0.7410];
  
elseif strcmp(colourString, 'azure')
        
    colour       = [0.3010, 0.7450, 0.9330];
    
elseif strcmp(colourString, 'yellow')
        
    colour       = [0.95,   0.85,   0.2];
    
elseif strcmp(colourString, 'orange')
        
    colour       = [1,      0.55,   0.1250];
    
elseif strcmp(colourString, 'violet')
        
    colour       = [0.4940, 0.1840, 0.5560];
    
elseif strcmp(colourString, 'granate')
        
    colour       = [0.6350, 0.0780, 0.1840]; 
    
elseif strcmp(colourString, 'gray')
        
    colour       = [0.6,    0.6,    0.6];
    
elseif strcmp(colourString, 'pink')
        
    colour       = [1,      0.4,    0.8];
    
elseif strcmp(colourString, 'peach')
        
    colour       = [1,      0.6,    0.6];
    
elseif strcmp(colourString, 'darkViolet')
        
    colour       = [0.3,    0,      0.6];

    
elseif strcmp(colourString, 'lightOrangeMap')
        
    colour       = [ones(1,256);     1:-(0.6)/255:0.4;   1:-1/255:0]';
    
elseif strcmp(colourString, 'lightBlueMap')
           
    colour       = [1:-0.5/255:0.5;  1:-0.8/255:0.2;     ones(1,256)]';

elseif strcmp(colourString, 'legendBackgroundColour')
    
     colour      = [1 1 1 0.75]';
    
end 
    

