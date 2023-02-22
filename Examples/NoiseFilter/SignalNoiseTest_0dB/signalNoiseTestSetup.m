% SETTINGS FOR THE MODE EXTRACTOR

% COPYRIGHT: 
%
%   MIT License
%   © 2023 Dario Chiantello <pipesoundlibrary@gmail.com>
%

% NOTES: Add this file in the [simulation folder] for the noise filter test
%
% OUTPUT:
%
%  - sns: signal to noise test setup structure
%

function sns = signalNoiseTestSetup()

sns = [];

% NAMES

sns.names.backgroundsFolder                         = 'Taps\Tap_1\Background_Observations';
sns.names.sourceNoiseRecording                      = 'Background_0';


% 
% %GENERAL
% 
% sns.general.sampleRate                              = 48000; 
% sns.general.recInitialTrimTime_s                    = 0.1;   
% 
% 
% 
% %SYNTHETIC SIGNAL 
% 
% sns.synth.duration_s                                = 5;                        %Default value in spaceSettings.m  
% sns.synth.timing_s                                  = [1,4];                    %Default value in spaceSettings.m      
% sns.synth.centreFrequencies_Hz                      = [1000 6000];              %Default value in spaceSettings.m  
% sns.synth.bandPassBands_Hz                          = [1000 2500];              %Default value in spaceSettings.m  
sns.synth.snRatio_dB                                = 0;                        %Default value in spaceSettings.m          
% sns.synth.randomSelection                           = false;                    %Default value in spaceSettings.m  
% 
% 
% 
% %NOISE MASK 
% 
% sns.filt.mask.disjointSet                           = true;                     %See getNoiseFrequencyMask.m, default value in spaceSettings.m
% sns.filt.mask.checkSet_pc                           = 25;                       %See getNoiseFrequencyMask.m, default value in spaceSettings.m
% sns.filt.mask.randomSelection                       = false;                    %See getNoiseFrequencyMask.m, default value in spaceSettings.m
% sns.filt.mask.windowDuration_s                      = 0.05;                     %See getNoiseFrequencyMask.m, default value in spaceSettings.m
% sns.filt.mask.windowOverlap_pc                      = 50;                       %See getNoiseFrequencyMask.m, default value in spaceSettings.m
% sns.filt.mask.noiseExcessThresholdSteps             = [1, 2, 5]; %[0:0.05:1]    %See getNoiseFrequencyMask.m, default value in spaceSettings.m
% sns.filt.mask.stdIncreaseStep                       = 0.05;                     %See getNoiseFrequencyMask.m, default value in spaceSettings.m
% 
% 
% 
% %FILTER 
% 
% sns.filt.filt.freqDensityLength_Hz                  = 200;                      %See noiseFilter.m, default value in spaceSettings.m          
% sns.filt.filt.timeDensityLength_s                   = 0.4;                      %See noiseFilter.m, default value in spaceSettings.m  
% sns.filt.filt.freqDensityThresholdMin               = 0.30;                     %See noiseFilter.m, default value in spaceSettings.m  
% sns.filt.filt.freqDensityThresholdMax               = 0.75;                     %See noiseFilter.m, default value in spaceSettings.m  
% sns.filt.filt.timeDensityThresholdMin               = 0.3;                      %See noiseFilter.m, default value in spaceSettings.m  
% sns.filt.filt.timeDensityThresholdMax               = 0.3;                      %See noiseFilter.m, default value in spaceSettings.m  
% sns.filt.filt.freqSmoothLengthSteps_Hz              = [300 500 1000];           %See noiseFilter.m, default value in spaceSettings.m                               	
% sns.filt.filt.timeSmoothLength_s                    = 0.25;                     %See noiseFilter.m, default value in spaceSettings.m  
% sns.filt.filt.firFilterOrder                        = 200;                      %See noiseFilter.m, default value in spaceSettings.m  
% 
% 
% 
% %FILTER OUTPUT
% 
% sns.filt.output.frameTime_s                         = 2.5;                      %See noiseFilter.m, default value in spaceSettings.m 
% sns.filt.output.firFilterOrder_low                  = 39;                       %See noiseFilter.m, default value in spaceSettings.m 
% sns.filt.output.firFilterOrder_high                 = 499;                      %See noiseFilter.m, default value in spaceSettings.m 
% sns.filt.output.noiseExcessThreshold                = 1;                        %See noiseFilter.m, default value in spaceSettings.m 
% sns.filt.output.freqSmoothLength_Hz                 = 300;                      %See noiseFilter.m, default value in spaceSettings.m 
% sns.filt.output.singleConv                          = false;                    %See noiseFilter.m, default value in spaceSettings.m 
% 
% 
% 
% %PERFORMANCE TABLES
% 
% sns.perfTable.freqSmoothLengthSteps_Hz              = [300, 500];               %Default value in spaceSettings.m 
% sns.perfTable.noiseExcessThresholdSteps             = [1, 2,  5];               %Default value in spaceSettings.m 
% 
% 
% 
% %COMPARISON TABLES
% 
% sns.compTable.freqSmoothLength_Hz                   = 500;                      %Default value in spaceSettings.m 
% sns.compTable.singleConv                            = 1;                        %Default value in spaceSettings.m 
% 
% 
% 
% %PICTURES
% 
% sns.pictures.format                                 = '';
% sns.pictures.size                                   = '';
% 


