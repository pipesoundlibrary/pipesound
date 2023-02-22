% DEFINITION OF NAMES AND PATHS

% COPYRIGHT: 
%
%   MIT License
%   © 2023 Dario Chiantello <pipesoundlibrary@gmail.com>
%


% PATHS

rootPath                                        = pwd;
codePath                                        = rootPath; 
utilitiesPath                                   = [rootPath,        '\', 'Utilities']; 
utilitiesCprintfPath                            = [utilitiesPath,   '\', 'Cprintf'];
acousticModelPath                               = [rootPath,        '\', 'AcousticModel']; 
noiseFilterPath                                 = [rootPath,        '\', 'NoiseFilter'];  



% SOUND REPOSITORY

ns.repository.path                              = 'E:\PipeSound\Examples';
        


% SOUNDBANK

ns.soundbank.folderName                         = 'SoundbankSnippet';               
ns.soundbank.backgrounds.observations.fileName  = 'background';    



% ACOUSTIC MODEL   

ns.acMod.folderName                             = 'AcousticModels';

ns.acMod.modes.folderName                       = 'Modes';
ns.acMod.modes.modeFolderRootName               = 'Mode';  
ns.acMod.modes.evanescentFolderName             = 'Evanescent'; 
ns.acMod.modes.propagativeFolderName            = 'Propagative';

ns.acMod.setup.folderName                       = 'SimulationSetup';
ns.acMod.setup.simulation.fileName              = 'simulationSetup';
ns.acMod.setup.kzSolver.fileName                = 'solverSetup';
ns.acMod.setup.extractor.fileName               = 'extractorSetup';

ns.acMod.kzSolver.resultRootName                = 'solverResult';
ns.acMod.extractor.cutFreqFileName              = 'cutFrequecies';
ns.acMod.extractor.resultsFileName              = 'modes';
ns.acMod.extractor.figureRootName               = 'modes';
ns.acMod.extractor.outputDataFileName           = 'nModesTable';



% FILTER AND SIGNAL TO NOISE RATIO TEST

ns.filter.folderName                            = 'NoiseFilter';

% filter
ns.filter.filterSummaryFileName                 = 'filterSummary'; 
ns.filter.filterModelFigureName                 = 'filterModel';    
ns.filter.sourceVsFilteredPicFileName           = 'sourceVsFiltered';   
ns.filter.windowPicFileName                     = 'window';
ns.filter.filteredSuffix                        = 'filtered';  
ns.filter.unfilteredSuffix                      = 'unfiltered'; 

% S/N test
ns.filter.test.settingsFileName                 = 'signalNoiseTestSetup';
ns.filter.test.syntheticSignalFigureName        = 'syntheticSignal'; 
ns.filter.test.sourceSignalsFigureName          = 'sourceSignals'; 
ns.filter.test.cleanSignalFileName              = 'cleanSignal';                  
ns.filter.test.noiseSignalFileName              = 'noiseSignal';    
ns.filter.test.noisySignalFileName              = 'noisySignal';    
ns.filter.test.cleanSignalFolderName            = 'FilteredCleanSignal';                    
ns.filter.test.noiseSignalFolderName            = 'FilteredNoiseSignal';                 
ns.filter.test.noisySignalFolderName            = 'FilteredNoisySignal';    
ns.filter.test.summaryFileName                  = 'summary';                     



