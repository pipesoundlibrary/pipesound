% SETTINGS FOR THE MODE EXTRACTOR

% COPYRIGHT: 
%
%   MIT License
%   © 2023 Dario Chiantello <pipesoundlibrary@gmail.com>
%

% NOTES: Add this file in [simulation folder]/SimulationSetup
%
% OUTPUT:
%
%  - exs: extractor settings structure
%

function exs = extractorSetup()

exs                                         = [];               % init struct in case of no fields


exs.settings.nMin                           = 0;                % See modeExtractor.m
exs.settings.nMax                           = 3;                % See modeExtractor.m
exs.settings.evanProp                       = 2;                % See modeExtractor.m
exs.settings.maxContinuityRadius_rad_m      = 1;                % See modeExtractor.m
exs.settings.differentialAverageSteps       = 10;               % See modeExtractor.m
exs.settings.minBranchLen                   = 8;                % See modeExtractor.m
exs.settings.initialGapMinK_rad_m           = 20;               % See modeExtractor.m
exs.settings.initialGapMinFreq_Hz           = 6000;             % See modeExtractor.m

%plot
exs.pictures.fLim_Hz                        = [];               % See modeExtractor.m
exs.pictures.kLim_rad_m                     = [];               % See modeExtractor.m
exs.pictures.colours                        = {};               % See modeExtractor.m
exs.pictures.format                         = '';               % See Utilities/getPlotSettings.m 
exs.pictures.size                           = '';               % See Utilities/getPlotSettings.m



%n=0
exs.pictures.colours.propagative{1,1}       = 'black';
exs.pictures.colours.propagative{1,2}       = 'violet';
exs.pictures.colours.propagative{1,3}       = 'orange';
exs.pictures.colours.propagative{1,4}       = 'blue';
exs.pictures.colours.propagative{1,5}       = 'lightGreen';
exs.pictures.colours.propagative{1,6}       = 'yellow';
exs.pictures.colours.propagative{1,7}       = 'azure';
exs.pictures.colours.propagative{1,8}       = 'pink';
exs.pictures.colours.propagative{1,9}       = '';
exs.pictures.colours.propagative{1,10}      = '';

exs.pictures.colours.evanescent{1,1}        = 'orange';
exs.pictures.colours.evanescent{1,2}        = 'blue';
exs.pictures.colours.evanescent{1,3}        = 'lightGreen';
exs.pictures.colours.evanescent{1,4}        = 'yellow';
exs.pictures.colours.evanescent{1,5}        = 'azure';
exs.pictures.colours.evanescent{1,6}        = 'pink';
exs.pictures.colours.evanescent{1,7}        = 'peach';
exs.pictures.colours.evanescent{1,8}        = 'darkViolet';
exs.pictures.colours.evanescent{1,9}        = '';
exs.pictures.colours.evanescent{1,10}       = '';


%n=1
exs.pictures.colours.propagative{2,1}       = 'red';
exs.pictures.colours.propagative{2,2}       = 'black';
exs.pictures.colours.propagative{2,3}       = 'violet';
exs.pictures.colours.propagative{2,4}       = 'orange';
exs.pictures.colours.propagative{2,5}       = 'blue';
exs.pictures.colours.propagative{2,6}       = 'lightGreen';
exs.pictures.colours.propagative{2,7}       = 'yellow';
exs.pictures.colours.propagative{2,8}       = 'azure';
exs.pictures.colours.propagative{2,9}       = 'pink';
exs.pictures.colours.propagative{2,10}      = '';

exs.pictures.colours.evanescent{2,1}        = 'granate';
exs.pictures.colours.evanescent{2,2}        = 'orange';
exs.pictures.colours.evanescent{2,3}        = 'blue';
exs.pictures.colours.evanescent{2,4}        = 'lightGreen';
exs.pictures.colours.evanescent{2,5}        = 'yellow';
exs.pictures.colours.evanescent{2,6}        = 'azure';
exs.pictures.colours.evanescent{2,7}        = 'pink';
exs.pictures.colours.evanescent{2,8}        = 'peach';
exs.pictures.colours.evanescent{2,9}        = 'darkViolet';
exs.pictures.colours.evanescent{2,10}       = '';


%n=2
exs.pictures.colours.propagative{3,1}       = 'red';
exs.pictures.colours.propagative{3,2}       = 'black';
exs.pictures.colours.propagative{3,3}       = 'violet';
exs.pictures.colours.propagative{3,4}       = 'orange';
exs.pictures.colours.propagative{3,5}       = 'blue';
exs.pictures.colours.propagative{3,6}       = 'lightGreen';
exs.pictures.colours.propagative{3,7}       = 'azure';
exs.pictures.colours.propagative{3,8}       = 'pink';
exs.pictures.colours.propagative{3,9}       = '';
exs.pictures.colours.propagative{3,10}      = '';

exs.pictures.colours.evanescent{3,1}        = 'black';
exs.pictures.colours.evanescent{3,2}        = 'yellow';
exs.pictures.colours.evanescent{3,3}        = 'violet';
exs.pictures.colours.evanescent{3,4}        = 'orange';
exs.pictures.colours.evanescent{3,5}        = 'blue';
exs.pictures.colours.evanescent{3,6}        = 'lightGreen';
exs.pictures.colours.evanescent{3,7}        = 'azure';
exs.pictures.colours.evanescent{3,8}        = 'pink';
exs.pictures.colours.evanescent{3,9}        = 'peach';
exs.pictures.colours.evanescent{3,10}       = 'darkViolet';
exs.pictures.colours.evanescent{3,11}       = 'gray';

%n=3
exs.pictures.colours.propagative{4,1}       = 'red';
exs.pictures.colours.propagative{4,2}       = 'black';
exs.pictures.colours.propagative{4,3}       = 'violet';
exs.pictures.colours.propagative{4,4}       = 'gray';
exs.pictures.colours.propagative{4,5}       = 'orange';
exs.pictures.colours.propagative{4,6}       = 'blue';
exs.pictures.colours.propagative{4,7}       = 'lightGreen';
exs.pictures.colours.propagative{4,8}       = 'azure';
exs.pictures.colours.propagative{4,9}       = 'pink';
exs.pictures.colours.propagative{4,10}      = '';

exs.pictures.colours.evanescent{4,1}        = 'black';
exs.pictures.colours.evanescent{4,2}        = 'yellow';
exs.pictures.colours.evanescent{4,3}        = 'granate';
exs.pictures.colours.evanescent{4,4}        = 'blue';
exs.pictures.colours.evanescent{4,5}        = 'lightGreen';
exs.pictures.colours.evanescent{4,6}        = 'azure';
exs.pictures.colours.evanescent{4,7}        = 'pink';
exs.pictures.colours.evanescent{4,8}        = 'peach';
exs.pictures.colours.evanescent{4,9}        = 'darkViolet';
exs.pictures.colours.evanescent{4,10}       = '';



