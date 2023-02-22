% SETTINGS FOR THE GEOMETRY AND THE MATERIALS OF THE SIMULATION

% COPYRIGHT: 
%
%   MIT License
%   © 2023 Dario Chiantello <pipesoundlibrary@gmail.com>
%

% NOTES: Add this file in [simulation folder]/SimulationSetup
%
% OUTPUT:
%
%  - sis: simulation settings structure
%

function sis = simulationSetup()

sis             = [];           % init struct in case of no fields;

%GEOMETRY
sis.R1_m        = 0.0984;       % See kzSolver.m
sis.R2_m        = 0.1000;       % See kzSolver.m


%MATERIAL PROPERTIES

%LIQUID - WATER
sis.c_phi_m_s   = 1479;         % See kzSolver.m
sis.ro_l_kg_m3  = 1000;         % See kzSolver.m

%SOLID SHIELD - ALUMINIUM
sis.c_gamma_m_s = 6420;         % See kzSolver.m  
sis.c_psi_m_s   = 3040;         % See kzSolver.m   
sis.ro_s_kg_m3  = 2700;         % See kzSolver.m


