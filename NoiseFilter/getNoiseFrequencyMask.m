% NOISE FILTER MASK

% COPYRIGHT: 
%
%   MIT License
%   © 2023 Dario Chiantello <pipesoundlibrary@gmail.com>
%

% DESCRIPTION:
% 
%   This script generates the noise mask used by noiseFilter.m. The noise 
%   mask is synthesised in the frequency domain by averaging the PSD of a
%   number of reference noise recordings in the characterisation-set 
%   and comparing this average against the equivalent obtained from the 
%   check-set. The noise mask is obtained as the level above the mean of 
%   the characterisation-set (mean + stdExcessFactor*std) that returns the 
%   maximum percentage (excessThreshold_pc) of check-set noise bins above 
%   the mask level. Background recordings must be organised as 
%   inputFileName.wav files contained in subfolders placed inside 
%   inputFolder. keepOutSubfolders are not considered for the creation 
%   of the mask. Any other content inside each subfolder is ignored. 
%

% INPUTS:
%
%   The input struct (maskArgs) is organised as described below. The 
%   necessary fields not specified are retrieved from spaceSettings.m as
%   default values. The fields inputFolder and inputFileName cannot be 
%   omitted.
%
%   maskArgs:
%     inputFolder           Path of the main folder cotaining the noise recordings. E.g.: 'E:\repository\backgrounds'
%     keepOutSubfolders     List of subfodlers to exclude. E.g.: {'background1', 'background2'} 
%     inputFileName         File name of the background recordings. E.g. 'background.wav'
%     disjointSet           True to keep characterization-set and check-set separate. False otherwise. Default:true 
%     checkSet_pc           Percentage of the background recordings to use as the check-set
%     randomSelection       True to divide randomly the set of background recordings. False otherwise. Default:true 
%     sampleRate            Audio sample rate. Audio are resampled using this value. Default see spaceSettings.m
%     recInitialTrimTime_s  Initial time length of the noise recordings to trim in seconds
%     windowDuration_s      Length of the windowing function in seconds. E.g.: 0.05
%     windowOverlap_pc      Overlap percentage of the windowing function. When not 50%, the window is modified to satisfy COLA. E.g.: 30, default: 50%.
%     excessThreshold_pc    Average percentace of check bins above the mask. E.g.: 1
%     stdIncreaseStep       Increment step of the noise mask STD. E.g. 0.05 
%     verboise              Display messages. true or false
%

% OUTPUTS:
%
%   noiseMask.
%      sampleRate           Sample rate of the aufio file     
%      mask.
%         mask              Noise mask obtained (mask = maskMean + stdMultiplier * maskStd)
%         mean              Noise mask mean   
%         std               Noise mask standard deviation
%         max               Noise mask max
%         stdExcessFactor   Std multiplier of the noise mask
%      win.
%         window            Source window
%         colaWindow        Modified COLA window used for the STFT
%         overlap           Number of overlapping samples between windows
%      outputString         Short summary of the mask
%

% PROCESSING STEPS: 
%
%   1- The background recordings are partitioned in a characterization-set 
%      and a check-set. The two sets can be the same (disjointSet = false) 
%      or disjoint. When the two partitions are disjoint, the selection can
%      also be performed randomly. The number of recordings in the two sets
%      is defined by checkSetPercentage.
% 
%   2- Slice the backgrounds in overlapping windows and extract the 
%      statistics from the one-sided PSD (mean, std, max).
% 
%   3- Calculate the mask mean and the mask std by averaging the mean and 
%      the std over the whole characterization-set. The mask max is the max
%      calculated over the whole characterization-set. The length of the 
%      mask is equal to the length of the one-sided PSD.
% 
%   4- Calculate the amplitude of the mask starting from the mean value and
%      increasing step by step by a fraction of std. The step is defined by 
%      stdIncreaseStep. The iteration stops when the average percentage 
%      of frequency bins of the check-set below the mask level is smaller 
%      than the excess threshold.
%

% NOTES:
%
%   The background recordings are resampled at the specified value. 
%   The same value of the sampleRate is also used in noiseFilter.m to 
%   resample the related signals.
% 

        
function noiseMask  = getNoiseFrequencyMask(maskArgs)                                
    
    %LOAD ENVIRONMENT VARIABLES

    nameSpace
    spaceSettings


    %SETTINGS AND DEFAULT VALUES

    if ~isfield(maskArgs, 'inputFolder') || ~isfield(maskArgs, 'inputFileName')
        error('Define a path for the background recordings and a name for the files')
    end
    
    if ~isfield(maskArgs, 'keepOutSubfolders')
        maskArgs.keepOutSubfolders = {};
    end  
     
    if ~isfield(maskArgs, 'sampleRate')
        maskArgs.sampleRate = ss.defaultSampleRate;
    end
    
    if ~isfield(maskArgs, 'recInitialTrimTime_s')
        maskArgs.recInitialTrimTime_s = ss.rec.recInitialTrimTime_s;
    end
    
    if ~isfield(maskArgs, 'verboise')
        maskArgs.verboise = false;
    end
    
    maskSettings.inputFolder            = [];
    maskSettings.inputFileName          = [];
    maskSettings.keepOutSubfolders      = [];
    maskSettings.disjointSet            = [];
    maskSettings.checkSet_pc            = [];
    maskSettings.randomSelection        = [];
    maskSettings.sampleRate             = [];
    maskSettings.recInitialTrimTime_s   = [];
    maskSettings.windowDuration_s       = [];
    maskSettings.windowOverlap_pc       = [];
    maskSettings.excessThreshold_pc     = [];
    maskSettings.stdIncreaseStep        = [];
    maskSettings.verboise               = [];
    
    if ~isfield(ss.filter, 'filterMask')
        ss.filt.filterMask = [];
    end
    
    maskSettings = getSettings(maskSettings, maskArgs, ss.filter.filterMask, false);
    
    if maskSettings.windowOverlap_pc > 50 || maskSettings.windowOverlap_pc < 25
        error('Set the window overlap between 25% and 50%');
    end 
    
    
    
    % BACKGROUNDS SELECTION AND PARTITION

    backgroundsFolder           = maskSettings.inputFolder;
    fileName                    = maskSettings.inputFileName;

    %get the background subfolders
    backgroundFoldersNames    	= getFolderNames(backgroundsFolder);

    removeItem = zeros(length(backgroundFoldersNames),1);

    for i = 1:length(maskSettings.keepOutSubfolders)

        itemIndex = findItemIndex(backgroundFoldersNames, maskSettings.keepOutSubfolders{i});

        removeItem(itemIndex) = 1;

    end

    backgroundFoldersNames = backgroundFoldersNames(~removeItem);

   

    numberOfBackgrounds = length(backgroundFoldersNames);

    %separate in characterization set and check set
    if maskSettings.disjointSet == false || numberOfBackgrounds < 4
    
        characterizationSet	= 1:numberOfBackgrounds;
        checkSet         	= characterizationSet;
        
    else
        
        if maskSettings.randomSelection

            randBackIndex      	= randperm(numberOfBackgrounds);
            characterizationSet	= randBackIndex(1:round(numberOfBackgrounds*(1-maskSettings.checkSet_pc/100)));    
            checkSet           	= randBackIndex(length(characterizationSet)+1 : end);

        else

            characterizationSet	= 1:round(numberOfBackgrounds*(1-maskSettings.checkSet_pc/100));
            checkSet         	= length(characterizationSet)+1 : numberOfBackgrounds;

        end
        
    end
    
    
    
    %GET COLA WINDOW
    
    windowLen           = roundEven(maskSettings.windowDuration_s * maskSettings.sampleRate);
    ovelappingSamples  	= floor(windowLen*maskSettings.windowOverlap_pc/100);
    [window, colaWin]   = colaWindow(windowLen, ovelappingSamples);
    
    

    %EXTRACT STATISTICS FROM BACKGROUNDS
    
    recordingsStatistics    = cell(length(backgroundFoldersNames),4);

    psdSamples              = floor(windowLen/2) + 1;
        
    recInitalTrimSamples    = maskSettings.recInitialTrimTime_s * maskSettings.sampleRate;
    
    for backgroundCounter = 1 : numberOfBackgrounds

        audioFolder = [backgroundsFolder, '\', backgroundFoldersNames{backgroundCounter}];

        %read audio files
        try
            background          = readAudio(audioFolder, fileName, maskSettings.sampleRate);
            background          = background(recInitalTrimSamples:end);
            
        catch
            error('.wav audio file missing in %s', audioFolder)
        end


        %obtaining fft transformed windows
        [winDoubleSidefft, winFreqBoundaries, winCentralTime] = ...
                        stft(   background,...
                                maskSettings.sampleRate,...
                                'Window',           colaWin,...
                                'OverlapLength',    ovelappingSamples,...
                                'FFTLength',        windowLen);

        %power spectral density            
        psdFull                 = 1/(maskSettings.sampleRate * windowLen) * abs(winDoubleSidefft(:,:)).^2;
        %take the second half including the Nyquist frequency
        psdHalf                 = psdFull(floor(windowLen/2):end,:);     
        %Pa^2/Hz - add the second half without freq 0 e nyquist freq  
        psdHalf(2:end-1,:)      = 2*psdHalf(2:end-1,:);        
        %dB/Hz Re 1Pa 
        psdHalfdB               = 10*log10(psdHalf);                             

        recordingsStatistics(backgroundCounter, 1) = {psdHalfdB};
        %db mean
        recordingsStatistics(backgroundCounter, 2) = {mean(psdHalfdB,2)};
        %db std
        recordingsStatistics(backgroundCounter, 3) = {std(psdHalfdB,0,2)};
        %db max
        recordingsStatistics(backgroundCounter, 4) = {max(psdHalfdB,[],2)};

    end


    %CREATE MASK

    %calculate average, std and max over the characterization set
    maskMean = zeros(psdSamples, 1);
    maskStd  = zeros(psdSamples, 1);
    maskMax  = zeros(psdSamples, 1) - 200;%init negative

    for backgroundCounter = 1 : length(characterizationSet)

        maskMean = maskMean +   recordingsStatistics{characterizationSet(backgroundCounter), 2};
        maskStd  = maskStd  +   recordingsStatistics{characterizationSet(backgroundCounter), 3};
        maskMax  = max(maskMax, recordingsStatistics{characterizationSet(backgroundCounter), 4});

    end

    maskMean = maskMean/length(characterizationSet);
    maskStd  = maskStd/length(characterizationSet);

    %add a proper margin to the mean mask using the check recordings
    mask                = maskMean;
    increaseCounter     = 0;

    %init to start the cycle
    averageExcessBins   = 100;

    while averageExcessBins > maskSettings.excessThreshold_pc

        averageExcessBins = 0;

        stdExcessFactor = maskSettings.stdIncreaseStep * increaseCounter;

        mask = maskMean + stdExcessFactor * maskStd;

        %check the number of std above the mean to include the vast majority of check noise
        for backgroundCounter = 1 : length(checkSet)

            checkNoise =  recordingsStatistics{checkSet(backgroundCounter), 1};

            checkMask  = repmat(mask,1,size(checkNoise,2)); %repeat per each window of the recorded noise

            checkTest  = checkNoise > checkMask;

            excessBins = sum(sum(checkTest))/(size(checkNoise,1)*size(checkNoise,2))*100;

            averageExcessBins = averageExcessBins + excessBins;
            
        end

        averageExcessBins = averageExcessBins / length(checkSet);
        increaseCounter   = increaseCounter + 1;

    end

    
    maskSummary = createMaskSummary(maskSettings, stdExcessFactor, averageExcessBins);
    
    if maskSettings.verboise == true

        fprintf(strrep(maskSummary,'\','\\'));
        
    end
    
                  
    noiseMask.sampleRate            = maskSettings.sampleRate;                 
    noiseMask.mask.mask             = mask;
    noiseMask.mask.mean             = maskMean;
    noiseMask.mask.std              = maskStd;
    noiseMask.mask.max              = maskMax;
    noiseMask.mask.stdExcessFactor  = stdExcessFactor;
    noiseMask.win.window            = window;
    noiseMask.win.colaWindow        = colaWin;
    noiseMask.win.overlap           = ovelappingSamples;
    noiseMask.outputString          = maskSummary;
 
end




%FUNCTIONS


% Check for change of sign given by discontinuities and remove false zeros in zeroCheckMatrix
%
% Inputs: 
%         zeroVector: bool col vector indication the position of the zeros;
%         modulusVector: col vector of abs(det(A));
%         discontinuityRadius: radius for the local maximum check
%
% Output: 
%         zeroVector: same as input but cleaned from false zeros

function [windowFunc, colaWin] = colaWindow(windowLen, ovelappingSamples)

    windowFunc                          = hann(windowLen);
    
    halfWindowLen                       = (round(windowLen/2));

    windowsSum                          = zeros(2*windowLen - ovelappingSamples, 1);
    windowsSum(1:windowLen)             = windowsSum(1:windowLen) + windowFunc;
    windowsSum(end-windowLen+1:end)     = windowsSum(end-windowLen+1:end) + windowFunc;

    diff                                = ones(halfWindowLen,1);
    diff(end-ovelappingSamples+1:end)   = (ovelappingSamples-1:-1:0)/(ovelappingSamples-1);

    unitaryMissing                      = 1 - windowsSum(round(windowLen/2)+1 : windowLen);

    modifiedWindow                      = zeros(size(windowFunc));
    modifiedWindow(halfWindowLen+1:end) = windowFunc(halfWindowLen+1:end) + diff(1:end).*unitaryMissing;
    modifiedWindow(1:halfWindowLen)     = modifiedWindow(end:-1:halfWindowLen+1);

    colaWin                             = modifiedWindow;

end



% Check for change of sign given by discontinuities and remove false zeros in zeroCheckMatrix
%
% Inputs: 
%         zeroVector: bool col vector indication the position of the zeros;
%         modulusVector: col vector of abs(det(A));
%         discontinuityRadius: radius for the local maximum check
%
% Output: 
%         zeroVector: same as input but cleaned from false zeros

function maskSummary = createMaskSummary(maskSettings, stdExcessFactor, averageExcessBins)

    inputFolderString       = sprintf("Backgrounds source folders = %s\n", maskSettings.inputFolder);
    keepOutString           = sprintf("Keepout backgrounds:\n");
    for i = 1:length(maskSettings.keepOutSubfolders)
        keepOutString       = strcat(keepOutString,  sprintf("   %s\n", maskSettings.keepOutSubfolders{i}));
    end
    inputFileString         = sprintf("Background file name: %s\n", maskSettings.inputFileName);
    randomSelectionString   = sprintf("Random selection: %d\n", maskSettings.randomSelection);
    disjointSetsString      = sprintf("Disjoint sets: %d\n", maskSettings.disjointSet);
    checkSetString          = sprintf("Check set: %d (percent)\n", maskSettings.checkSet_pc);
    initialTrimString       = sprintf("Initial trim time: %.3f (s)\n", maskSettings.recInitialTrimTime_s);
    windowDurationString    = sprintf("Window duration: %.3f (s)\n", maskSettings.windowDuration_s);
    windowOverlapString     = sprintf("Window overlap: %.1f (percent)\n", maskSettings.windowOverlap_pc);
    stdIncreaseString       = sprintf("Std increase step: %.3f\n", maskSettings.stdIncreaseStep);
    stdString               = sprintf("STD above the mean: %.2f \n", stdExcessFactor);
    noiseExcessString       = sprintf("Noise excess: %.3f (percent)\n", averageExcessBins);

    maskSummary = strcat(   inputFolderString,...
                            keepOutString,...
                            inputFileString,...
                            randomSelectionString,...
                            disjointSetsString,...
                            checkSetString,...
                            initialTrimString,...
                            windowDurationString,...
                            windowOverlapString,...
                            stdIncreaseString,...
                            stdString,...
                            noiseExcessString );

end

