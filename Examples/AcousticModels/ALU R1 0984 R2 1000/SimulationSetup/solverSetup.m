% SETTINGS FOR THE KZ SOLVER

% COPYRIGHT: 
%
%   MIT License
%   © 2023 Dario Chiantello <pipesoundlibrary@gmail.com>
%

% NOTES: Add this file in [simulation folder]/SimulationSetup
%
% INPUTS (only to change settings in runtime):  
%       
% - n: n mode order
% - evanescent: true for evanescent modes, false for propagative modes
%
% OUTPUT:
%
%  - sos: solver settings structure
%

function sos = solverSetup(n, evanescent)

sos                                             = [];       % init struct in case of no fields

sos.settings.nMin                               = 0;        % see kzSolver.m
sos.settings.nMax                               = 0;        % see kzSolver.m

sos.settings.evanProp                           = 2;        % see kzSolver.m

sos.settings.kMin_rad_m                         = 0;        % see kzSolver.m              
sos.settings.kMax_rad_m                         = 75;       % see kzSolver.m
sos.settings.kResolution_rad_m                  = 0.01;     % see kzSolver.m   
sos.settings.kSearchRadius_rad_m                = 1.5;      % see kzSolver.m
sos.settings.kDiscontinuityRadius_rad_m         = 0.08;     % see kzSolver.m

sos.settings.cutFrequencySearchScale            = 10;       % see kzSolver.m
sos.settings.fDiscontinuityRadius_Hz            = 10;       % see kzSolver.m

sos.settings.fMin_Hz                            = 5;        % see kzSolver.m           
sos.settings.fMax_Hz                            = 25000;    % see kzSolver.m            
sos.settings.fResolution_Hz                     = 5;        % see kzSolver.m 
sos.settings.fInitSteps                         = 15;       % see kzSolver.m
sos.settings.fZeroNeighboursSteps               = 4;        % see kzSolver.m
sos.settings.zeroCutExtraInitFreqLeap_Hz        = 50;       % see kzSolver.m
sos.settings.zeroCutExtraInitPointsRadius       = 4;        % see kzSolver.m
sos.settings.fExtraInitSteps_Hz                 = 25000;    % see kzSolver.m (e.g. [1035, 1498])

sos.pictures.format                             = '';       % See Utilities/getPlotSettings.m
sos.pictures.size                               = '';       % See Utilities/getPlotSettings.m


% if evanescent
% 
%     sos.settings.kMax_rad_m                   = 500; 
% 
% else
% 
%     sos.settings.kMax_rad_m                   = 250; 
% 
% end