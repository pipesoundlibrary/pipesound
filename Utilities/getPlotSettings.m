% GET PLOT SETTINGS

% COPYRIGHT: 
%
%   MIT License
%   © 2023 Dario Chiantello <pipesoundlibrary@gmail.com>
%

% DESCRIPTION: 
%
%   Set format, size and properties of the output pictures. If the struct
%   userSettings is null, returns the default values.
%

% INPUTS:
%   userSettings.
%      format                   CUSTOM/IEEE
%      size                     DEFAULT/FULL_1/THREEQUARTERS_1/HALF_1

% OUTPUTS:
%   plotSettings.
%      paperSize
%      position
%      
%      lines.
%         thin
%         thick
%
%      labels.
%         fontSize
%         tickFontSize
%
%      legend.
%         fontSize
%         boxLineThickness
%         boxColourType
%         boxColourData
%


function plotSettings = getPlotSettings(userSettings)    

    if ~exist('userSettings', 'var')
        userSettings.format = ''; 
        userSettings.size = ''; 
    end

    if ~isfield(userSettings, 'format') 
        userSettings.format = ''; 
    end
        
    if ~isfield(userSettings, 'size') 
        userSettings.size = ''; 
    end
        
    plotSettings.lines.thin                 = getPlotLineThickness('thinLine');
    plotSettings.lines.thick                = getPlotLineThickness('thickLine');
    
    plotSettings.legend.boxLineThickness    = getPlotLineThickness('defaultPlotLineThickness');
    plotSettings.legend.boxColourType       = 'truecoloralpha';
    plotSettings.legend.boxColourData       = uint8(255*getColour('legendBackgroundColour')); 
    
    
    %IEEE STYLE
	if strcmp(userSettings.format, 'IEEE')
    
        plotSettings.labels.fontSize            = getPlotFontSize('ieeeLabelFontSize');
        plotSettings.labels.tickFontSize        = getPlotFontSize('defaultTickFontsize');
        
        plotSettings.legend.fontSize            = getPlotFontSize('ieeeTickFontSize');
 
        if     strcmp(userSettings.size, 'FULL_1')
        
            plotSettings.position               = getPlotDimensions('181x565_figurePosition');
            plotSettings.paperSize              = getPlotDimensions('181x56_paperSize');
        
        elseif strcmp(userSettings.size, 'THREEQUARTERS_1')
        
            plotSettings.position               = getPlotDimensions('160x56_figurePosition');
            plotSettings.paperSize              = getPlotDimensions('160x56_paperSize');
        
        elseif strcmp(userSettings.size, 'HALF_1')
            
            plotSettings.position               = getPlotDimensions('88x50_figurePosition');
            plotSettings.paperSize              = getPlotDimensions('88x50_paperSize');
        
        else %default
            
            plotSettings.position               = getPlotDimensions('160x56_figurePosition');
            plotSettings.paperSize              = getPlotDimensions('160x56_paperSize');
            
        end
        
                
    %CUSTOM STYLE   
	else
                  
        plotSettings.labels.fontSize            = getPlotFontSize('defaultLabelFontsize');
        plotSettings.labels.tickFontSize        = getPlotFontSize('defaultTickFontsize');
        
        plotSettings.legend.fontSize            = getPlotFontSize('defaultLegendFontsize');
        
        
        if     strcmp(userSettings.size, 'FULL_1')
        
            plotSettings.position               = getPlotDimensions('165x52_figurePosition');
            plotSettings.paperSize              = getPlotDimensions('165x52_paperSize');
        
        elseif strcmp(userSettings.size, 'THREEQUARTERS_1')
            
            plotSettings.position               = getPlotDimensions('125x67_figurePosition');
            plotSettings.paperSize              = getPlotDimensions('125x67_paperSize');
        
        elseif strcmp(userSettings.size, 'HALF_1')
        
            plotSettings.position               = getPlotDimensions('65x45_figurePosition');
            plotSettings.paperSize              = getPlotDimensions('65x45_paperSize');

        else %default
        
            plotSettings.labels.fontSize        = getPlotFontSize('default');
            %plotSettings.labels.tickFontSize    = getPlotFontSize('default');
            plotSettings.legend.fontSize        = getPlotFontSize('default');
            
            plotSettings.position               = getPlotDimensions('default_figurePosition');
            plotSettings.paperSize              = getPlotDimensions('default_paperSize');
            
        end
            
	end
       
end





