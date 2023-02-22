% GENERAL SETTINGS AND DEFAULT VALUES

% COPYRIGHT: 
%
%   MIT License
%   © 2023 Dario Chiantello <pipesoundlibrary@gmail.com>
%


% GENERAL DEFAULT SETTINGS
ss.defaultSampleRate                                    = roundEven(48000);                             %default sampling rate value (samples/second)
ss.defaultWindowDuration_s                              = 0.05;



% RECORDINGS
ss.rec.recInitialTrimTime_s                             = 0.1;                                          



% ACOUSTIC MODEL

%simulation settings
ss.acMod.simSettings.R1_m                               = 0.0492;                                       %default value. See kzSolver.m
ss.acMod.simSettings.R2_m                               = 0.0508;                                       %default value. See kzSolver.m
ss.acMod.simSettings.c_phi_m_s                          = 1479;                                         %default value. See kzSolver.m 
ss.acMod.simSettings.ro_l_kg_m3                         = 1000;                                         %default value. See kzSolver.m
ss.acMod.simSettings.c_gamma_m_s                        = 6420;                                         %default value. See kzSolver.m   
ss.acMod.simSettings.c_psi_m_s                          = 3040;                                         %default value. See kzSolver.m 
ss.acMod.simSettings.ro_s_kg_m3                         = 2700;                                         %default value. See kzSolver.m

%solver
ss.acMod.solver.nMin                                    = 0;                                            %default value. See kzSolver.m
ss.acMod.solver.nMax                                    = 0;                                            %default value. See kzSolver.m
ss.acMod.solver.evanProp                                = 2;                                            %default value. See kzSolver.m
ss.acMod.solver.kMin_rad_m                              = 0;                                            %default value. See kzSolver.m
ss.acMod.solver.kMax_rad_m                              = 75;                                           %default value. See kzSolver.m
ss.acMod.solver.kResolution_rad_m                       = 0.01;                                         %default value. See kzSolver.m
ss.acMod.solver.kSearchRadius_rad_m                     = 1.5;                                          %default value. See kzSolver.m
ss.acMod.solver.kDiscontinuityRadius_rad_m              = 0.08;                                         %default value. See kzSolver.m
ss.acMod.solver.cutFrequencySearchScale                 = 10;                                           %default value. See kzSolver.m
ss.acMod.solver.fDiscontinuityRadius_Hz                 = 10;                                           %default value. See kzSolver.m
ss.acMod.solver.fMin_Hz                                 = 5;                                            %default value. See kzSolver.m
ss.acMod.solver.fMax_Hz                                 = 25000;                                        %default value. See kzSolver.m
ss.acMod.solver.fResolution_Hz                          = 5;                                            %default value. See kzSolver.m
ss.acMod.solver.fInitSteps                              = 15;                                           %default value. See kzSolver.m
ss.acMod.solver.fZeroNeighboursSteps                    = 4;                                            %default value. See kzSolver.m
ss.acMod.solver.zeroCutExtraInitFreqLeap_Hz             = 50;                                           %default value. See kzSolver.m
ss.acMod.solver.zeroCutExtraInitPointsRadius            = 4;                                            %default value. See kzSolver.m
ss.acMod.solver.fExtraInitSteps_Hz                      = 25000;                                        %default value. See kzSolver.m

%extractor
ss.acMod.extractor.nMin                                 = 0;                                            %default value. See modeExtractr.m
ss.acMod.extractor.nMax                                 = 0;                                            %default value. See modeExtractr.m
ss.acMod.extractor.evanProp                             = 2;                                            %default value. See modeExtractr.m
ss.acMod.extractor.maxContinuityRadius_rad_m            = 1;                                            %default value. See modeExtractr.m
ss.acMod.extractor.differentialAverageSteps             = 10;                                           %default value. See modeExtractr.m
ss.acMod.extractor.minBranchLen                         = 8;                                            %default value. See modeExtractr.m
ss.acMod.extractor.initialGapMinK_rad_m                 = 20;                                           %default value. See modeExtractr.m
ss.acMod.extractor.initialGapMinFreq_Hz                 = 6000;                                         %default value. See modeExtractr.m



% NOISE FILTER DEFAULT VALUES

%noise mask
ss.filter.filterMask.disjointSet                        = true;                                         %default value. See getNoiseFrequencyMask.m
ss.filter.filterMask.checkSet_pc                        = 25;                                           %default value. See getNoiseFrequencyMask.m
ss.filter.filterMask.randomSelection                    = true;                                         %default value. See getNoiseFrequencyMask.m
ss.filter.filterMask.windowDuration_s                   = ss.defaultWindowDuration_s;                   %default value. See getNoiseFrequencyMask.m
ss.filter.filterMask.windowOverlap_pc                   = 50;                                           %default value. See getNoiseFrequencyMask.m
ss.filter.filterMask.excessThreshold_pc                 = 1;                                            %default value. See getNoiseFrequencyMask.m
ss.filter.filterMask.stdIncreaseStep                    = 0.05;                                         %default value. See getNoiseFrequencyMask.m
 
%filter
ss.filter.filter.freqDensityLength_Hz                   = 200;                                          %default value. See noiseFilter.m
ss.filter.filter.timeDensityLength_s                    = 0.4;                                          %default value. See noiseFilter.m
ss.filter.filter.freqDensityThresholdMin                = 0.30;                                         %default value. See noiseFilter.m
ss.filter.filter.freqDensityThresholdMax                = 0.75;                                         %default value. See noiseFilter.m
ss.filter.filter.timeDensityThresholdMin                = 0.30;                                         %default value. See noiseFilter.m
ss.filter.filter.timeDensityThresholdMax                = 0.30;                                         %default value. See noiseFilter.m
ss.filter.filter.freqSmoothLength_Hz                    = 300;                                          %default value. See noiseFilter.m
ss.filter.filter.timeSmoothLength_s                     = 0.25;                                         %default value. See noiseFilter.m
ss.filter.filter.firFilterOrder                         = 200;                                          %default value. See noiseFilter.m
ss.filter.filter.firFilterOrder_low                     = 40;                                           %default value. See noiseFilter.m
ss.filter.filter.firFilterOrder_high                    = 250;                                          %default value. See noiseFilter.m
ss.filter.filter.singleLinearConvolution                = false;                                        %default value. See noiseFilter.m

%signal-to-noise test - synthetic signal
ss.filter.snTest.synth.duration_s                       = 5;                                            %default value. See signalNoiseRatioTestUtility.m
ss.filter.snTest.synth.timing_s                         = [1,4];                                        %default value. See signalNoiseRatioTestUtility.m
ss.filter.snTest.synth.centreFrequencies_Hz             = [1000 6000];                                  %default value. See signalNoiseRatioTestUtility.m
ss.filter.snTest.synth.bandPassBands_Hz                 = [1000 2500];                                  %default value. See signalNoiseRatioTestUtility.m
ss.filter.snTest.synth.snRatio_dB                       = 0;                                            %default value. See signalNoiseRatioTestUtility.m
ss.filter.snTest.synth.randomSelection                  = false;                                        %default value. See signalNoiseRatioTestUtility.m

%signal-to-noise test - noise mask settings
ss.filter.snTest.mask.disjointSet                       = ss.filter.filterMask.disjointSet;             %default value. See getNoiseFrequencyMask.m, signalNoiseRatioTestUtility.m
ss.filter.snTest.mask.checkSet_pc                       = ss.filter.filterMask.checkSet_pc;             %default value. See getNoiseFrequencyMask.m, signalNoiseRatioTestUtility.m
ss.filter.snTest.mask.randomSelection                   = false;                                        %default value. See getNoiseFrequencyMask.m, signalNoiseRatioTestUtility.m
ss.filter.snTest.mask.windowDuration_s                  = ss.filter.filterMask.windowDuration_s;        %default value. See getNoiseFrequencyMask.m, signalNoiseRatioTestUtility.m
ss.filter.snTest.mask.windowOverlap_pc                  = ss.filter.filterMask.windowOverlap_pc;        %default value. See getNoiseFrequencyMask.m, signalNoiseRatioTestUtility.m
ss.filter.snTest.mask.noiseExcessThresholdSteps         = [1, 2, 5];                                    %default value. See getNoiseFrequencyMask.m, signalNoiseRatioTestUtility.m
ss.filter.snTest.mask.stdIncreaseStep                   = ss.filter.filterMask.stdIncreaseStep;         %default value. See getNoiseFrequencyMask.m, signalNoiseRatioTestUtility.m

%signal-to-noise test - filter settings
ss.filter.snTest.filt.freqDensityLength_Hz              = ss.filter.filter.freqDensityLength_Hz;        %default value. See noiseFilter.m, signalNoiseRatioTestUtility.m
ss.filter.snTest.filt.timeDensityLength_s               = ss.filter.filter.timeDensityLength_s;         %default value. See noiseFilter.m, signalNoiseRatioTestUtility.m
ss.filter.snTest.filt.freqDensityThresholdMin           = ss.filter.filter.freqDensityThresholdMin;     %default value. See noiseFilter.m, signalNoiseRatioTestUtility.m
ss.filter.snTest.filt.freqDensityThresholdMax           = ss.filter.filter.freqDensityThresholdMax;     %default value. See noiseFilter.m, signalNoiseRatioTestUtility.m
ss.filter.snTest.filt.timeDensityThresholdMin           = ss.filter.filter.timeDensityThresholdMin;     %default value. See noiseFilter.m, signalNoiseRatioTestUtility.m
ss.filter.snTest.filt.timeDensityThresholdMax           = ss.filter.filter.timeDensityThresholdMax;     %default value. See noiseFilter.m, signalNoiseRatioTestUtility.m
ss.filter.snTest.filt.freqSmoothLengthSteps_Hz          = [300 500 1000];                               %default value. See noiseFilter.m, signalNoiseRatioTestUtility.m
ss.filter.snTest.filt.timeSmoothLength_s                = ss.filter.filter.timeSmoothLength_s;          %default value. See noiseFilter.m, signalNoiseRatioTestUtility.m
ss.filter.snTest.filt.firFilterOrder                    = ss.filter.filter.firFilterOrder;              %default value. See noiseFilter.m, signalNoiseRatioTestUtility.m

%signal-to-noise test - filter output
ss.filter.snTest.output.frameTime_s                     = 2.5;                                          %default value. See noiseFilter.m, signalNoiseRatioTestUtility.m
ss.filter.snTest.output.firFilterOrder_low              = 39;                                           %default value. See noiseFilter.m, signalNoiseRatioTestUtility.m
ss.filter.snTest.output.firFilterOrder_high             = 499;                                          %default value. See noiseFilter.m, signalNoiseRatioTestUtility.m
ss.filter.snTest.output.noiseExcessThreshold            = 1;                                            %default value. See noiseFilter.m, signalNoiseRatioTestUtility.m
ss.filter.snTest.output.freqSmoothLength_Hz             = 300;                                          %default value. See noiseFilter.m, signalNoiseRatioTestUtility.m
ss.filter.snTest.output.singleConv                      = false;                                        %default value. See noiseFilter.m, signalNoiseRatioTestUtility.m

%signal-to-noise test - output performance table
ss.filter.snTest.perfTable.freqSmoothLengthSteps_Hz     = [300, 500];                                   %default value. See signalNoiseRatioTestUtility.m
ss.filter.snTest.perfTable.noiseExcessThresholdSteps    = [1, 2,  5];                                   %default value. See signalNoiseRatioTestUtility.m

%signal-to-noise test - output comparison table
ss.filter.snTest.compTable.freqSmoothLength_Hz          = 500;                                          %default value. See signalNoiseRatioTestUtility.m
ss.filter.snTest.compTable.singleConv                   = 1;                                            %default value. See signalNoiseRatioTestUtility.m





