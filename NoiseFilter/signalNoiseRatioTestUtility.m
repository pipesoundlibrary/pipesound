% UTILITY FOR FILTER PERFORMANCE ASSESSMENT USING A SYNTHETIC SIGNAL

% COPYRIGHT: 
%
%   MIT License
%   © 2023 Dario Chiantello <pipesoundlibrary@gmail.com>
%

% DESCRIPTION:
% 
%   This utility is designed to test the performance of the filter by
%   measurin, signal attenuation, noise attenuation, and signal-to-noise
%   ratio for different configurations of the filter. The test signal 
%   is an artificial signal composed of real background noise and an 
%   pure artificial signal obtained by filtering white noise in 
%   user-defined frequency ranges. Background noise and pure signals are 
%   processed separately using the same filter synthesised using the noisy 
%   artificial signal as the synthesis reference. Processing background 
%   noise and pure signals separately allows precise calculation of 
%   filter attenuations and S/N before and after the filter. 
%

% INPUTS:
%
%   The setting for a specific test should be defined in a configuration 
%   file signalNoiseTestSetup.m stored in the test folder. Examples of 
%   the configuration file can be found in Examples\NiseFilter. Fields not
%   specified in the configuration file (except for the sns.names struct) 
%   are retrieved from spaceSettings.m as default values. The data 
%   structures of the configuration file are listed below, along with the 
%   related description.
%
%
%   signalNoiseTestSetup.m
%
%   sns.names.                      Paths and names
%      backgroundsFolder            Path of the soundbank containing the background noise recordings.
%      sourceNoiseRecording         Folder name of the noise used for the synthetic signal (not for the filter).
%
%   sns.general.                    General settings
%     sampleRate                    Signals sample rate (samples/sec) 
%     recInitialTrimTime_s          Time to skip at the beginning of the noise recordings  
%
%   sns.synth.                      Synthetic test signal
%     duration_s                    Test signal duration in seconds.
%     timing_s                      Start and stop time of the clean signal in seconds. E.g.: [1 4]
%     centreFrequencies_Hz          Centre frequencies of the band-pass filters for the synthesis of the clean signal in Hz. E.g.: [2000 4000]             
%     bandPassBands_Hz              Bands of the band-pass filters for the synthesis of clean signal in Hz. E.g.: [1000 500] 
%     snRatio_dB                    Desired signal to noise ratio of the synthetic signal in dB. E.g.: -20
%     randomSelection               Select a random chunck from the background. True or false. If false, select the first part. 
%
%   sns.simSettings.                Filter settings
%     noiseExcessThresholdSteps     Steps of noise excess threshold. See getNoiseFilterMask.m. E.g.: [0:0.05:1, 1.1:0.1:1.9, 2:0.5:5].
%     disjointSet                   See getNoiseFilterMask.m
%     maskRandomSelection           See getNoiseFilterMask.m
%     firFilterOrder                Order of the filter used. See noiseFilter.m
%     freqDensityLength_Hz          See noiseFilter.m
%     timeDensityLength_s           See noiseFilter.m
%     freqDensityThresholdMin       See noiseFilter.m   
%     freqDensityThresholdMax       See noiseFilter.m  
%     timeDensityThresholdMin       See noiseFilter.m  
%     timeDensityThresholdMax       See noiseFilter.m
%     timeSmoothLength_s            See noiseFilter.m
%     freqSmoothLengthSteps_Hz      Steps of the freqSmoothLength_Hz. See noiseFilter.m. E.g.: [300 500 1000]; 
%
%   sns.filterExampleOutput.        Settings used only for the filter selected for the output results 
%     frameTime_s                   Central time of the frame in seconds used to generate the filter model pictures. E.g.: 2.5.
%     firFilterOrder_low            Low-order filter used for filter model illustration. E.g.: 100.
%     firFilterOrder_high           High-order filter used for filter model illustration. E.g.: 500.
%     noiseExcessThreshold          Selected noise excess threshold percentage of the noise mask. E.g.:1. See getNoiseFilterMask.m
%     freqSmoothLength_Hz           Selected length of the smooth filter. See noiseFilter.m 
%     singleConv                    True for forward convolution only. False for forward/backward convolution. See noiseFilter.m.
%
%   sns.perfTableOutput.            Setting used to generate the performance table
%     noiseExcessThresholdSteps     Noise excess threshold steps of the noise mask. See getNoiseFilterMask.m. E.g.:[0.5, 1, 3].
%     freqSmoothLengthSteps_Hz      Frequecny smooth length steps. See noiseFilter.m. E.g.:[50, 200].
%
%   sns.compTableOutput.            Settings used to generate the comparison table
%     freqSmoothLength_Hz           Frequency smooth length. See noiseFilter.m. E.g.:50.
%     singleConv                    True for forward convolution only. False for forward/backward convolution only. See noiseFilter.m.
%
%   sns.pictures.                   Output picture settings. 
%     format = '';                  Picture format. See getPlotSettings.m
%     size   = '';                  Picture size. See getPlotSettings.m
%

% OUTPUTS:
%
%   Synthetic signal:
%
%     cleanSignal.wav               Synthetic clean signal audio file.
%     noiseSignal.wav               Background noise audio file.
%     noisySignal.wav               Synthetic signal audio file.
%     sourceSignals.pdf             Noise and clean signal (time domain) used to create the synthetic signal.
%     syntheticSignal.pdf           Time domain and PSD of the synthetic signal.       
%
%   Filter (for the filter and the frame selected):
%  
%     Filter output for clean signal (audio, filter model figure, source vs fitered figure, filterSummary.txt)
%     Filter output for noise signal (audio, filter model figure, source vs fitered figure, filterSummary.txt)
%     Filter output for noisy signal (audio, filter model figure, source vs fitered figure, filterSummary.txt)
%
%   summary.txt                     Filter and performance summary
%

% PROCESSING STEPS: 
%
%   1- Creation of a synthetic test signal. Three signals are generated: 
%      cleanSignal, noiseSignal, and noisySignal. The latter is simply the 
%      superimposition of the former two. The clean signal is created by
%      band-pass filtering white noise. The noiseSignal is extracted from
%      the background recordings. The noise recording selected is excluded
%      from the synthesis of the noise mask in step 2. The desired
%      signal-to-noise ratio is imposed by changing the amplitude of the
%      clean signal with respect to the amplitude of the background. The
%      energy is calculated over the whole duration of the signals. All the
%      generated signals are saved as .wav for reference and comparison
%      with different software libraries. 
%  
%   2- Noise filtering. Noise filtering is repeated for multiple values of
%      noise excess threshold, smoothing lengths, and for forward and
%      forward/backward convolution. The noise mask is synthesised using
%      the background recordings (except the one used for the synthesis of 
%      the artificial signal). The filter synthesis reference is the 
%      synthetic noisySignal. All the signals (cleanSignal, noiseSignal, 
%      noisySignal) are processed using the same filter in order to measure
%      the new signal-to-noise ratio precisely. A set of results is 
%      generated for the desired filter setting and for the desired frame.
%      A summary text file reports the settings and the results of the 
%      tests.
%
%
% NOTES:
%
%   The clean signal is generated by filtering white noise. Therefore,
%   a different version of the clean signal is obtained every time the
%   utility is launched, even using the same parameters. Results slightly
%   different should be expected when repeating the same test.
%   To compare the results against other filter solutions, use the copy of
%   the signals saved in the simulation folder.
%


clear
clc


%DEFAULT SIMULATION NAME
defaultTestFolderName  = 'SignalNoiseTest_0dB';


%LOAD LIBRARY

loadLibray();
nameSpace
spaceSettings



%NAME OF THE SIMULATION FOLDER

prompt          = sprintf("Default test folder name: ""%s"". Press enter to confirm or type a new name.\n", defaultTestFolderName);
testFolderName  = input(prompt, "s");

if isempty(testFolderName)
    
    testFolderName = defaultTestFolderName;
    
end



%SAVE/OVERWRITE DATA

prompt         = sprintf("Previous result data in the same folder will be overwritten. Do you want to continue [Y/N] >\n");
overwriteAns   = input(prompt, "s");
if ~strcmpi(overwriteAns, 'y')
    return
end



%SETTINGS - FOLDERS AND FILES NAMES


%test folder
names.signalNoiseTestFolder = [ns.repository.path, '\', ns.filter.folderName, '\', testFolderName];
if ~exist(names.signalNoiseTestFolder, 'dir')
   mkdir(names.signalNoiseTestFolder)
end

%simulation setup
settingsFile = [names.signalNoiseTestFolder,'\',ns.filter.test.settingsFileName, '.m'];
if ~exist(settingsFile, 'file')
    error('Add a setup file in the folder: %s.', names.signalNoiseTestFolder)     
end
addpath(names.signalNoiseTestFolder);
sns = signalNoiseTestSetup();

if isfield(sns, 'names')
    
    %background folder
    if isfield(sns.names, 'backgroundsFolder')
        names.backgroundsFolder = [ns.repository.path, '\',  ns.soundbank.folderName, '\', sns.names.backgroundsFolder];
    else
        error('names.backgroundsFolder not defined in %s', settingsFile)
    end
    
    if isempty(sns.names.backgroundsFolder) || ~exist(names.backgroundsFolder, 'dir')
        error('names.backgroundsFolder defined in %s is pointing to a missing folder', settingsFile)
    end
    
    
    %sourceNoiseRecording
    if isfield(sns.names, 'sourceNoiseRecording')
        names.sourceNoiseRecording = sns.names.sourceNoiseRecording;
    else
        error('names.sourceNoiseRecording not defined in %s', settingsFile)
    end
    
    if isempty(sns.names.sourceNoiseRecording)
        error('Invalid definition of names.sourceNoiseRecording in %s', settingsFile)
    end
    
    
else
    error('names struct not defined in %s', settingsFile)
end
    
names.noiseFileName = [ns.soundbank.backgrounds.observations.fileName, '.wav'];



%SETTINGS - GENERAL

general.sampleRate              = ss.defaultSampleRate;
general.recInitialTrimTime_s    = ss.rec.recInitialTrimTime_s;

if isfield(sns, 'general')

    if isfield(sns.general, 'sampleRate')
        general.sampleRate = sns.general.sampleRate;
    end
    
    if isfield(sns.general, 'recInitialTrimTime_s')
        general.recInitialTrimTime_s = sns.general.recInitialTrimTime_s;
    end
    
end


%SETTINGS - SYNTHETIC SIGNAL

synth.duration_s                            = [];
synth.timing_s                              = [];
synth.centreFrequencies_Hz                  = []; 
synth.bandPassBands_Hz                      = [];
synth.snRatio_dB                            = [];
synth.randomSelection                       = [];

if ~isfield(sns, 'synth')
    sns.synth = [];
end
if ~isfield(ss.filter.snTest, 'synth')
    ss.filter.snTest = [];
end
synth = getSettings(synth, sns.synth, ss.filter.snTest.synth, false);



%SETTINGS - NOISE MASK

simSettings.mask.disjointSet                = [];
simSettings.mask.checkSet_pc             	= [];
simSettings.mask.randomSelection            = [];
simSettings.mask.windowDuration_s           = [];
simSettings.mask.windowOverlap_pc           = [];
simSettings.mask.noiseExcessThresholdSteps  = [];
simSettings.mask.stdIncreaseStep            = [];

if ~isfield(sns, 'filt')
    sns.filt = [];
end
if ~isfield(sns.filt, 'mask')
    sns.filt.mask = [];
end
if ~isfield(ss.filter.snTest, 'mask')
    ss.filter.snTest.mask = [];
end
simSettings.mask = getSettings(simSettings.mask, sns.filt.mask, ss.filter.snTest.mask, false);



%SETTINGS - NOISE FITLER

simSettings.filt.freqDensityLength_Hz       = [];            
simSettings.filt.timeDensityLength_s        = []; 
simSettings.filt.freqDensityThresholdMin    = [];     
simSettings.filt.freqDensityThresholdMax    = [];     
simSettings.filt.timeDensityThresholdMin    = [];          
simSettings.filt.timeDensityThresholdMax    = [];            
simSettings.filt.freqSmoothLengthSteps_Hz   = [];                                    	
simSettings.filt.timeSmoothLength_s         = [];         
simSettings.filt.firFilterOrder             = [];  

if ~isfield(sns.filt, 'filt')
    sns.filt.filt = [];
end
if ~isfield(ss.filter.snTest, 'filt')
    ss.filter.snTest.filt = [];
end
simSettings.filt = getSettings(simSettings.filt, sns.filt.filt, ss.filter.snTest.filt, false);



%SETTINGS - OUTPUT

simSettings.output.frameTime_s              = [];
simSettings.output.firFilterOrder_low       = [];
simSettings.output.firFilterOrder_high      = [];
simSettings.output.noiseExcessThreshold     = [];
simSettings.output.freqSmoothLength_Hz      = [];
simSettings.output.singleConv               = [];

if ~isfield(sns.filt, 'output')
    sns.filt.output = [];
end
if ~isfield(ss.filter.snTest, 'output')
    ss.filter.snTest.output = [];
end
simSettings.output = getSettings(simSettings.output, sns.filt.output, ss.filter.snTest.output, false);



%SETTINGS - PERFORMANCE TABLE

perfTableOutput.freqSmoothLengthSteps_Hz    = [];  
perfTableOutput.noiseExcessThresholdSteps   = [];

if ~isfield(sns, 'perfTable')
    sns.perfTable = [];
end
if ~isfield(ss.filter.snTest, 'perfTable')
    ss.filter.snTest.perfTable = [];
end
perfTableOutput = getSettings(perfTableOutput, sns.perfTable, ss.filter.snTest.perfTable, false);



%SETTINGS - COMPARISON TABLE

compTableOutput.freqSmoothLength_Hz         = [];
compTableOutput.singleConv                  = [];

if ~isfield(sns, 'compTable')
    sns.compTable = [];
end
if ~isfield(ss.filter.snTest, 'compTable')
    ss.filter.snTest.compTable = [];
end
compTableOutput = getSettings(compTableOutput, sns.compTable, ss.filter.snTest.compTable, false);



%SETTINGS - PICTURES

pictures.format = '';
pictures.size   = '';

if isfield(sns, 'pictures')

    if isfield(sns.pictures, 'format')
        pictures.format = sns.pictures.format;
    end
    
    if isfield(sns.pictures, 'size')
        pictures.size = sns.pictures.size;
    end    
    
end





%GET NOISE

sourceNoiseFolder = [names.backgroundsFolder, '\', names.sourceNoiseRecording];

%read noise audio file
try
    [noiseSignal, physicalUnit] = readAudio(sourceNoiseFolder, names.noiseFileName, general.sampleRate);
catch
    error('Audio file %s missing in %s', names.noiseFileName, sourceNoiseFolder)
end

%extract a random chunk of the desired length from the source noise 
numberOfSamples = synth.duration_s * general.sampleRate;
initTrimSamples = general.recInitialTrimTime_s * general.sampleRate;

if length(noiseSignal) - initTrimSamples < numberOfSamples
    error('The noise signal selected is too short')
else  
    
    if synth.randomSelection 
    
        startSample = randi([initTrimSamples+1, length(noiseSignal) - synth.duration_s*general.sampleRate + 1]);
        noiseSignal = noiseSignal(startSample: startSample + numberOfSamples -1); 

    else
        
        noiseSignal = noiseSignal(initTrimSamples+1: initTrimSamples + numberOfSamples); 
        
    end
end



%CREATE CLEAN SIGNAL

%check filter bands
if length(synth.centreFrequencies_Hz) ~= length(synth.bandPassBands_Hz) 
    error('Number of centre frequencies and bands must be the same')
end
    
bandPassMinFreq = synth.centreFrequencies_Hz - synth.bandPassBands_Hz/2;
bandPassMaxFreq = synth.centreFrequencies_Hz + synth.bandPassBands_Hz/2;

if min(bandPassMinFreq) < 0 ||  max(bandPassMaxFreq) > general.sampleRate/2 
    error('Required frequencies of the clean signal beyond the allowed spectrum')
end
       
%create bandpass filter
[filterPoints, filterShape]     = findFilterPoints(synth, general.sampleRate);
b                               = fir2(200, filterPoints, filterShape);

%create white noise and apply the band-pass filter
cleanSignal                     = rand(numberOfSamples,1)-0.5;
cleanSignal                     = real(filtfilt( b, 1, cleanSignal));

%limit the duration of the clean signal in time
timeWindowLen                   = round(numberOfSamples * (synth.timing_s(2) - synth.timing_s(1)) / synth.duration_s);
startSilenceLen                 = round((numberOfSamples - timeWindowLen)/2);
endSilenceLen                   = numberOfSamples - timeWindowLen - startSilenceLen;
    
cleanSignalWindow               = [zeros(1,startSilenceLen) hann(timeWindowLen)' zeros(1,endSilenceLen)];
cleanSignal                     = cleanSignal.*cleanSignalWindow';

%determine the amplitude to match the desired S/N
noiseEnergyDb                   = 10*log10(sum(noiseSignal.^2));

alpha                           = 10^((synth.snRatio_dB + noiseEnergyDb)/10);
beta                            = sum(cleanSignal.^2);

cleanSignal                     = sqrt(alpha/beta) * cleanSignal;
cleanSignalEnergyDb             = 10*log10(sum(cleanSignal.^2));


%CREATE NOISY SIGNAL

noisySignal                     = cleanSignal + noiseSignal;

%play noisy signal
soundsc(noisySignal, general.sampleRate)



%SAVE SYNTHETIC SIGNALS AND PLOTS

createSyntheticSignalPlots(pictures, noiseSignal, cleanSignal, noisySignal, general.sampleRate, physicalUnit, names.signalNoiseTestFolder)

noiseSignalFileName = [ns.filter.test.noiseSignalFileName, '_', ns.filter.unfilteredSuffix, '.wav'];
cleanSignalFileName = [ns.filter.test.cleanSignalFileName, '_', ns.filter.unfilteredSuffix, '.wav'];
noisySignalFileName = [ns.filter.test.noisySignalFileName, '_', ns.filter.unfilteredSuffix, '.wav'];

saveAudio(names.signalNoiseTestFolder, noiseSignalFileName, noiseSignal, general.sampleRate, physicalUnit, '');      
saveAudio(names.signalNoiseTestFolder, cleanSignalFileName, cleanSignal, general.sampleRate, physicalUnit, '');          
saveAudio(names.signalNoiseTestFolder, noisySignalFileName, noisySignal, general.sampleRate, physicalUnit, '');


% ASSESS FILTER PERFORMANCE

%noise mask static parameters
maskArgs.inputFolder                        = names.backgroundsFolder;
maskArgs.keepOutSubfolders                  = {names.sourceNoiseRecording};
maskArgs.inputFileName                      = names.noiseFileName; 
maskArgs.disjointSet                        = simSettings.mask.disjointSet;
maskArgs.checkSet_pc                        = simSettings.mask.checkSet_pc;
maskArgs.randomSelection                    = simSettings.mask.randomSelection;
maskArgs.sampleRate                         = general.sampleRate;
maskArgs.recInitialTrimTime_s               = general.recInitialTrimTime_s;
maskArgs.windowDuration_s                   = simSettings.mask.windowDuration_s;
maskArgs.windowOverlap_pc                   = simSettings.mask.windowOverlap_pc;
maskArgs.stdIncreaseStep                    = simSettings.mask.stdIncreaseStep; 
maskArgs.verboise                           = false; 

%filter static parameters
filterArgs.signals.synthesisReference	    = [names.signalNoiseTestFolder, '\' noisySignalFileName];

filterArgs.pars.freqDensityLength_Hz        = simSettings.filt.freqDensityLength_Hz;
filterArgs.pars.timeDensityLength_s         = simSettings.filt.timeDensityLength_s;
filterArgs.pars.freqDensityThresholdMin     = simSettings.filt.freqDensityThresholdMin; 
filterArgs.pars.freqDensityThresholdMax     = simSettings.filt.freqDensityThresholdMax; 
filterArgs.pars.timeDensityThresholdMin     = simSettings.filt.timeDensityThresholdMin; 
filterArgs.pars.timeDensityThresholdMax     = simSettings.filt.timeDensityThresholdMax;        
filterArgs.pars.timeSmoothLength_s          = simSettings.filt.timeSmoothLength_s;
filterArgs.pars.firFilterOrder              = simSettings.filt.firFilterOrder;

filterArgs.pictures.format                  = pictures.format;
filterArgs.pictures.size                    = pictures.size;
filterArgs.pictures.frameTime_s             = simSettings.output.frameTime_s;
filterArgs.pars.firFilterOrder_low          = simSettings.output.firFilterOrder_low;
filterArgs.pars.firFilterOrder_high         = simSettings.output.firFilterOrder_high;


noiseSignalFolder                           = [names.signalNoiseTestFolder, '\', ns.filter.test.noiseSignalFolderName];
cleanSignalFolder                           = [names.signalNoiseTestFolder, '\', ns.filter.test.cleanSignalFolderName];
noisySignalFolder                           = [names.signalNoiseTestFolder, '\', ns.filter.test.noisySignalFolderName];

%find the step where to create the output
[~, noiseExcessThresholdOutputStep]         = min(abs(simSettings.mask.noiseExcessThresholdSteps - simSettings.output.noiseExcessThreshold));
[~, freqSmoothOutputStep]                   = min(abs(simSettings.filt.freqSmoothLengthSteps_Hz  - simSettings.output.freqSmoothLength_Hz));
singleConvStep                              = simSettings.output.singleConv ;

performanceLabels   =  {    'noiseExcess',...
                            'stdExcessFactor',...
                            'singleConv',...
                            'freqSmooth',...
                            'snGain',...
                            'signalAttenuation',...
                            'noiseAttenuation'};

tableLength         = length(simSettings.mask.noiseExcessThresholdSteps) * length(simSettings.filt.freqSmoothLengthSteps_Hz) * 2;

performance         = table(zeros(tableLength,1),...
                            zeros(tableLength,1),...
                            zeros(tableLength,1),...
                            zeros(tableLength,1),...
                            zeros(tableLength,1),...
                            zeros(tableLength,1),...
                            zeros(tableLength,1),...
                            'VariableNames',performanceLabels);    
                        
                       
itemCounter = 0;

for excessCounter = 1:length(simSettings.mask.noiseExcessThresholdSteps)

    for singleConvCounter = 0:1
    
        for freqSmoothCounter = 1:length(simSettings.filt.freqSmoothLengthSteps_Hz)
    
            itemCounter = itemCounter + 1;

            %get noise mask
            maskArgs.excessThreshold_pc             = simSettings.mask.noiseExcessThresholdSteps(excessCounter);     
            noiseMask                               = getNoiseFrequencyMask(maskArgs);

            %filter signals
            filterArgs.signals.noiseMask           	= noiseMask;
            filterArgs.pars.singleLinearConvolution = singleConvCounter;
            filterArgs.pars.freqSmoothLength_Hz     = simSettings.filt.freqSmoothLengthSteps_Hz(freqSmoothCounter);
            
            %output filter results
            filterArgs.pictures.create              = excessCounter     == noiseExcessThresholdOutputStep &&...
                                                      singleConvCounter == singleConvStep &&...
                                                      freqSmoothCounter == freqSmoothOutputStep;
                                                         
            %filter signals
            filterArgs.signals.sourceSignal         = [names.signalNoiseTestFolder, '\' noiseSignalFileName];
            filterArgs.pictures.outputFolder        = noiseSignalFolder;
            filteredNoiseSignal                     = noiseFilter(filterArgs);
            filterArgs.signals.sourceSignal         = [names.signalNoiseTestFolder, '\' cleanSignalFileName];
            filterArgs.pictures.outputFolder        = cleanSignalFolder;
            filteredCleanSignal                     = noiseFilter(filterArgs);

            %save filtered signal
            if filterArgs.pictures.create
           
                filterArgs.signals.sourceSignal     = [names.signalNoiseTestFolder, '\'  noisySignalFileName];
                filterArgs.pictures.outputFolder    = [names.signalNoiseTestFolder, '\', ns.filter.test.noisySignalFolderName];
                [filteredNoisySignal, ~, ~, savedFilterString]...
                                                    = noiseFilter(filterArgs);                
                
                filteredNoiseSignalFileName         = [ns.filter.test.noiseSignalFileName, '_', ns.filter.filteredSuffix, '.wav'];
                filteredCleanSignalFileName       	= [ns.filter.test.cleanSignalFileName, '_', ns.filter.filteredSuffix, '.wav'];
                filteredNoisySignalFileName      	= [ns.filter.test.noisySignalFileName, '_', ns.filter.filteredSuffix, '.wav'];

                saveAudio(noiseSignalFolder, filteredNoiseSignalFileName, filteredNoiseSignal, general.sampleRate, physicalUnit, '');
                saveAudio(cleanSignalFolder, filteredCleanSignalFileName, filteredCleanSignal, general.sampleRate, physicalUnit, '');
                saveAudio(noisySignalFolder, filteredNoisySignalFileName, filteredNoisySignal, general.sampleRate, physicalUnit, '');

            end

            %calculate results
            filteredNoiseEnergyDb       = 10*log10(sum(filteredNoiseSignal.^2));
            filteredCleanSignalEnergyDb = 10*log10(sum(filteredCleanSignal.^2)); 

            snGain                      = filteredCleanSignalEnergyDb - filteredNoiseEnergyDb - synth.snRatio_dB;
            signalAttenuation           = cleanSignalEnergyDb - filteredCleanSignalEnergyDb;
            noiseAttenuation            = noiseEnergyDb       - filteredNoiseEnergyDb;

            
            %save results
            performance.noiseExcess(itemCounter)        = simSettings.mask.noiseExcessThresholdSteps(excessCounter);
            performance.stdExcessFactor(itemCounter)    = noiseMask.mask.stdExcessFactor;
            performance.singleConv(itemCounter)         = singleConvCounter;
            performance.freqSmooth(itemCounter)         = simSettings.filt.freqSmoothLengthSteps_Hz(freqSmoothCounter);
            performance.snGain(itemCounter)             = snGain;
            performance.signalAttenuation(itemCounter)  = signalAttenuation;
            performance.noiseAttenuation(itemCounter)   = noiseAttenuation;

        end
            
    end
        
end
    

%CREATE PERFORMANCE OUTPUT AND SUMMARY

%create performance table
performanceTable = createPerformanceTable(performance, simSettings, perfTableOutput);

%createComparisonTable
performanceComparisonTable = createComparisonTable(performance, simSettings, compTableOutput);

%create summary
createSummary(names.backgroundsFolder, names.sourceNoiseRecording, synth, savedFilterString, performanceTable, performanceComparisonTable, names.signalNoiseTestFolder);

fprintf("Done!\n")






%FUNCTIONS


% Load the libray
%

function loadLibray()

    rootPath            = mfilename('fullpath');
    indx                = strfind(rootPath, '\');
    indx                = indx(end-1);
    libraryRootPath     = rootPath(1:indx-1);
    cd(libraryRootPath);
    path(path, libraryRootPath);
    pathManager

end



% SYNTHESIS OF THE BANDPASS FILTER
%
%              *---------*
%             |           |
%            |             |
%   *-------*               *---------*
%  fp1     fp2 fp3     fp4 fp5       fp6

function [filterPoints, filterShape] = findFilterPoints(synth, sampleRate)

    spaceSettings

    nf                          = sampleRate/2;
    filterLen                   = synth.duration_s*nf + 1;
    
    filterPoints                = linspace(0,1,filterLen);
    filterShape                 = zeros(1, filterLen); 
    
    for filterCounter = 1:length(synth.centreFrequencies_Hz)
    
        f2 = synth.centreFrequencies_Hz(filterCounter) - synth.bandPassBands_Hz(filterCounter)/2;
        f5 = synth.centreFrequencies_Hz(filterCounter) + synth.bandPassBands_Hz(filterCounter)/2;
        f3 = f2 +   (f5 - f2)/5;
        f4 = f2 + 4*(f5 - f2)/5;
 
        f2p = round(f2/nf * filterLen);
        f3p = round(f3/nf * filterLen);
        f4p = round(f4/nf * filterLen);
        f5p = round(f5/nf * filterLen);
        
        bandFilter = zeros(1, filterLen); 
        
        bandFilter(f2p:f3p) = ((f2p:f3p) - f2p)/(f3p - f2p);
        bandFilter(f3p:f4p) = 1;
        bandFilter(f4p:f5p) = 1 - ((f4p:f5p) - f4p)/(f5p - f4p);
            
        filterShape = max(filterShape, bandFilter);
        
    end
    
end



%PLOT SYNTHETIC SIGNAL

function createSyntheticSignalPlots(pictures, noiseSignal, cleanSignal, noisySignal, sampleRate, unit, outputFolder)

    nameSpace

    settings.format                 = pictures.format;
    settings.size                   = pictures.size;
    plotSettings                    = getPlotSettings(settings);

    numberOfSamples                 = length(noisySignal);
    testSignalsDuration             = numberOfSamples/sampleRate;
    t                               = linspace(0,testSignalsDuration, numberOfSamples);
    f                               = linspace(0,sampleRate/2,floor(numberOfSamples/2)+1);

    
    %plot the noisy synthetic signal    
    noisySignalfft                  = fft(noisySignal);
    noisySignalPsd                  = 1/(sampleRate * length(noisySignalfft)) * abs(noisySignalfft).^2;
    noisySignalPsdHalf              = noisySignalPsd(1 : length(noisySignalfft)/2 +1);
    noisySignalPsdHalf(2:end-1)     = 2*noisySignalPsdHalf(2:end-1);
    noisySignalPsdHalfDb            = 10*log10(noisySignalPsdHalf);

    figureSyntheticSignal = figure('name', 'noisy signals');

    subplot(2,1,1)
    plot(t,noisySignal, 'LineWidth', plotSettings.lines.thin)
    
    xlim([0 t(end)])
    try
    ylim([-1.1 * max(abs(noisySignal))  1.1 * max(abs(noisySignal))])
    catch
    ylim([-1 1]) 
    end
    xlabel('$Time~[s]$',...
                'fontsize', plotSettings.labels.fontSize, 'interpreter','latex');
    ylabel(['$Amplitude~', unit, '$'],...
                'fontsize', plotSettings.labels.fontSize, 'interpreter','latex');
    set(gca,'fontsize', plotSettings.labels.tickFontSize);
            
            
            
    subplot(2,1,2)
    plot(f/1000, noisySignalPsdHalfDb, 'LineWidth', plotSettings.lines.thin)
    
    xlim([0 f(end)/1000])
    xlabel('$Frequency~[kHz]$',...
                'fontsize', plotSettings.labels.fontSize, 'interpreter','latex');
    ylabel(['$PSD~[dB/Hz~Re1', unit(2:end-1), ']$'],...
                'fontsize', plotSettings.labels.fontSize, 'interpreter','latex');
    set(gca, 'fontsize', plotSettings.labels.tickFontSize);  
            
            
    set(    figureSyntheticSignal,...
            'Units',        'centimeters',...
            'Position',     plotSettings.position,...
            'PaperUnits',   'centimeters',...
            'PaperSize',    plotSettings.paperSize);


    figureFile  = [outputFolder, '\', ns.filter.test.syntheticSignalFigureName, '.pdf']; 
    saveas(figureSyntheticSignal, figureFile);
        
    close(figureSyntheticSignal)

    

    %plot noise and clean signal
    figureSourceSignals = figure('name', 'signals: noise, clean');
    
    subplot(2,1,1)
    plot(t,noiseSignal)
    
    xlim([0 t(end)])
    try
    ylim([-1.1 * max(abs(noiseSignal))  1.1 * max(abs(noiseSignal))])
    catch
    ylim([-1 1]); 
    end
    
    xlabel('$Time~[s]$',...
                'fontsize', plotSettings.labels.fontSize, 'interpreter','latex');

    ylabel(['$Amplitude~', unit, '$'],...
                'fontsize', plotSettings.labels.fontSize, 'interpreter','latex');
    
    set(gca, 'fontsize', plotSettings.labels.tickFontSize);  
            
            
    lgdItemsN{1}            = 'Noise'; 

    [lgd,objh]              = legend(lgdItemsN, 'FontSize', plotSettings.legend.fontSize);    
    lgd.Location            = 'NorthEast';
    lgd.NumColumns          = 1;
    lgd.LineWidth           = plotSettings.legend.boxLineThickness;
    lgd.BoxFace.ColorType   = plotSettings.legend.boxColourType;
    lgd.BoxFace.ColorData   = plotSettings.legend.boxColourData;
            
    
    subplot(2,1,2)
    plot(t,cleanSignal)
    
    xlim([0 t(end)])
    try
    ylim([-1.1 * max(abs(cleanSignal))  1.1 * max(abs(cleanSignal))])
    catch
    ylim([-1 1]); 
    end
    
    xlabel('$Time~[s]$',...
                'fontsize', plotSettings.labels.fontSize, 'interpreter','latex');

    ylabel(['$Amplitude~', unit, '$'],...
                'fontsize', plotSettings.labels.fontSize, 'interpreter','latex');
    
    set(gca,'fontsize', plotSettings.labels.tickFontSize);
            
    lgdItemsS{1}            = 'Signal'; 

    [lgd,objh]              = legend(lgdItemsS, 'FontSize', plotSettings.legend.fontSize);   
    lgd.Location            = 'NorthEast';
    lgd.NumColumns          = 1;
    lgd.LineWidth           = plotSettings.legend.boxLineThickness;
    lgd.BoxFace.ColorType   = plotSettings.legend.boxColourType;
    lgd.BoxFace.ColorData   = plotSettings.legend.boxColourData;
            
       
    set(    figureSourceSignals,...
            'Units',        'centimeters',...
            'Position',     plotSettings.position,...
            'PaperUnits',   'centimeters',...
            'PaperSize',    plotSettings.paperSize);


    figureFile  = [outputFolder, '\', ns.filter.test.sourceSignalsFigureName, '.pdf']; 
    saveas(figureSourceSignals, figureFile);

    close(figureSourceSignals)

end



%PERFORMANCE TABLE

function performanceTable = createPerformanceTable(performance, simSettings, perfTableOutput)

    %performance table          
    performanceTable = array2table(zeros(0,size(performance, 2)), 'VariableNames', performance.Properties.VariableNames);

    for noiseExcessStep = 1:length(perfTableOutput.noiseExcessThresholdSteps)

        desiredNoiseExcessValue = perfTableOutput.noiseExcessThresholdSteps(noiseExcessStep);

        [~, noiseExcessValueId] = min(abs(simSettings.mask.noiseExcessThresholdSteps - desiredNoiseExcessValue));
        noiseExcessValue        = simSettings.mask.noiseExcessThresholdSteps(noiseExcessValueId);

        for freqSmoothStep = 1:length(perfTableOutput.freqSmoothLengthSteps_Hz)

            desiredFreqSmoothValue = perfTableOutput.freqSmoothLengthSteps_Hz(freqSmoothStep);

            [~, freqSmoothValueId] = min(abs(simSettings.filt.freqSmoothLengthSteps_Hz - desiredFreqSmoothValue));
            freqSmoothValue        = simSettings.filt.freqSmoothLengthSteps_Hz(freqSmoothValueId);

            performanceTable = [performanceTable;...
                                performance(...
                                performance.noiseExcess == noiseExcessValue &...
                                performance.freqSmooth == freqSmoothValue, :)];
        end
    end    
end

    
  
%COMPARISON TABLE

function performanceComparisonTable = createComparisonTable(performance, simSettings, compTableOutput)

    %comparison table                          
    desiredSmoothValue = compTableOutput.freqSmoothLength_Hz;
    [~, smoothValueId] = min(abs(simSettings.filt.freqSmoothLengthSteps_Hz - desiredSmoothValue));
    freqSmoothValue    = simSettings.filt.freqSmoothLengthSteps_Hz(smoothValueId);

    performanceComparisonTable = performance( performance.freqSmooth == freqSmoothValue &...
                                              performance.singleConv == compTableOutput.singleConv, :);
                                          
end



%SUMMARY

function createSummary(backgroundsFolder, sourceNoiseRecording, synth, filterString, performanceTable, performanceComparisonTable, signalNoiseTestFolder)

    nameSpace

    summaryFile = [signalNoiseTestFolder, '\', ns.filter.test.summaryFileName, '.txt'];
    
    fid=fopen(summaryFile,'w');

    printFileAndPrompt(fid, 'NOISE SOURCE\n', 'red')
    printFileAndPrompt(fid, ['Backgrounds folder: ', strrep(backgroundsFolder, '\', '\\'), '\n'])
    printFileAndPrompt(fid,'\n\n')
    
    
    
    printFileAndPrompt(fid, 'SYNTHETIC SIGNAL\n', 'red')
    printFileAndPrompt(fid, sprintf("Background noise file: %s\n", sourceNoiseRecording))
    printFileAndPrompt(fid, sprintf("Background random selection: %d\n", synth.randomSelection))
    printFileAndPrompt(fid, sprintf("Duration: %.3f (s)\n", synth.duration_s))
    printFileAndPrompt(fid, sprintf("Clean signal time: [%.3f, %.3f] (s)\n", synth.timing_s(1), synth.timing_s(2)))
    printFileAndPrompt(fid, sprintf("Clean signal components (Hz): \n"))
    
    for i=1:length(synth.centreFrequencies_Hz)
    
        printFileAndPrompt(fid, sprintf("   Centre: %.1f, Band: %.1f\n", synth.centreFrequencies_Hz(i), synth.bandPassBands_Hz(i)))
 
    end
    
    printFileAndPrompt(fid, sprintf("S/N dB: %.1f (dB)\n", synth.snRatio_dB))
    printFileAndPrompt(fid,'\n\n')
    
    printFileAndPrompt(fid, 'FILTER OF THE OUTPUT SAVED\n', 'red')
    printFileAndPrompt(fid, strrep(filterString, '\', '\\'))
    printFileAndPrompt(fid,'\n\n')
    
    printFileAndPrompt(fid, 'PERFORMANCE\n', 'red')
    
    performanceTableString              = writeTableString(performanceTable);
    performanceComparisonTableString    = writeTableString(performanceComparisonTable);  
    
    printFileAndPrompt(fid,'\n')
    printFileAndPrompt(fid, 'Performance table:\n', 'blue')
    printFileAndPrompt(fid, performanceTableString);
    
    printFileAndPrompt(fid,'\n')
    printFileAndPrompt(fid, 'Comparison table:\n', 'blue')
    printFileAndPrompt(fid, performanceComparisonTableString);

    fclose(fid);
    
end

%print table
function tableString = writeTableString(inputTable)

    tableString = [];

    stringItemLen = max(cellfun(@(x) length(x), inputTable.Properties.VariableNames)) + 3;
    
    for lineCounter = 1:size(inputTable,2)
        
        stringItem(1:stringItemLen) = ' ';
        
        varName = inputTable.Properties.VariableNames{lineCounter};
        stringItem(1:length(varName)) = varName;
    
        tableString = [tableString, sprintf(stringItem)];
        
    end
    
    tableString = [tableString, sprintf('\n')];    
    
    for lineCounter = 1:size(inputTable,1)
        for colCounter = 1:size(inputTable,2)
            
            stringItem(1:stringItemLen) = ' ';
        
            if colCounter == 3
                
                if inputTable{lineCounter, colCounter} == 1
                    
                    value = ' F';
                    
                else
                    
                    value = ' F/B';
                    
                end
                
            else
                
                value = sprintf(' %0.3f', inputTable{lineCounter, colCounter});
                
            end
            
            stringItem(1:length(value)) = value;

            tableString = [tableString, stringItem];
            
        end
        
        tableString = [tableString, sprintf('\n')];
        
    end

end






