% NOISE FILTER

% COPYRIGHT: 
%
%   MIT License
%   © 2023 Dario Chiantello <pipesoundlibrary@gmail.com>
%

% DESCRIPTION:
% 
%   The noise filter aims to attenuate the background noise of a source 
%   signal by synthesising a FIR filter per each processed window. The 
%   filter synthesis requires a synthesis reference signal and a noise mask 
%   obtained from getNoiseFrequencyMask.m. The window filters are obtained
%   by:
% 
%   - comparing the reference signal and the noise mask,
%   - applying frequency and time domain bin clustering, 
%   - applying a final frequency and time gaussian smoothing.
% 
%   Filtering is executed in the time domain with forward single linear or 
%   forward/backward double linear convolution. Filtered windows are then 
%   overlapped, accounting for the leading and trailing tails to minimise 
%   filtering artefacts. See filterExample.m for a simple example using the
%   SoundbankSnippet. 
%

% INPUTS: 
%
%   The input struct (filterArgs) is organised as described below. 
%   Necessary Fields not specified are retrieved from spaceSettings.m as 
%   default values. The fields noiseMask and sourceSignal in 
%   filterArgs.signals cannot be omitted. If the field synthesisReference
%   is not specified, the sourceSignal is used as the synthesisReference. 
%   See getNoiseFrequencyMask.m for the synthesis of the noiseMask. 
%
%   filterArgs.signals.
%     noiseMask                     Noise mask obtained from getNoiseFrequencyMask.m.
%     sourceSignal                  Full file path of the signal to filter. E.g: 'E:\datasets\test.wav'.
%     synthesisReference            Full file path of the synthesis reference signal (.wav or .mp4). Leave null unless running S/N test. Default: [].
%   filterArgs.pars.
%     freqDensityLength_Hz          Length of the frequency neighbourhood in Hz. E.g.: 600. 
%     timeDensityLength_s           Length of the time neighbourhood in seconds. E.g.: 0.5.
%     freqDensityThresholdMin       Min threshold of frequency neighbours. Range: 0:1, E.g.: 0.15. 
%     freqDensityThresholdMax       Max threshold of frequency neighbours. Range: 0:1, E.g.: 0.35. 
%     timeDensityThresholdMin       Min threshold of time neighbours. Range: 0:1, E.g.: 0.15.  
%     timeDensityThresholdMax       Max threshold of time neighbours. Range: 0:1, E.g.: 0.20. 
%     freqSmoothLength_Hz           Length in Hz of the gaussian window to smooth the mask in frequency. E.g. 300. 
%     timeSmoothLength_s            Length in seconds of the gaussian window to smooth the mask in time. E.g. 0.25.
%     firFilterOrder                Order of the synthesised FIR filter. E.g.: 50. 
%     firFilterOrder_low            Low-order filter used for filter model illustration. E.g.: 15.
%     firFilterOrder_high           High-order filter used for filter model illustration. E.g.: 200.
%     singleLinearConvolution       True for sigle linear convolution, false for double linear convolution. 
%   filterArgs.pictures.
%     create                        True to create the pictures, false otherwise.
%     frameTime_s                   Time frame of the filter illustrated in seconds. E.g: 3.45, deafult: half signal length.
%     format                        Format of the pictures. See getPlotSettings.m.
%     size                          Size of the pictures. See getPlotSettings.m. 
%     visible                       Set true to show the figures while being generated, false otherwise.  
%     plotWindow                    Set true to plot the windowing function. Deafult: false.
%     outputFolder                  Full path of the output folder to save audio and pictures. E.g: 'E:\testFolder', default: []
%   filterArgs.summary.
%     verboise                      True to display summary, false otherwise. 
%

% OUTPUTS:
%
%   variables:
%
%      filteredSignal               Filtered signal
%      sourceSignal                 Source signal
%      signalUnit                   Physical unit of the signal (from sourceSignal metadata)
%      outputString                 Short summary of the filter
%
%   files:
%
%      filterModelFigure.pdf        Illustration of the filter model.
%      sourceVsFilteredFigure.pdf   Plots comparing source and filtered signals.
%      windowFunctionFigure.pdf     Representation of the window function.
%      filterSummary.txt            Text file reporting the filter summary.
%

% PROCESSING STEPS: 
%
%   1- The signal to filter (sourceSignal) is loaded along with the noise 
%      mask (noiseMask) obtained from getNoiseFrequencyMask. The loaded 
%      signals are sliced into overlapping windows and time-smoothed using 
%      the windowing function (and the related overlap) imported from
%      getNoiseFrequencyMask. Windows and overlap should always satisfy the
%      COLA condition.
% 
%   2- The power spectral density of the windowed synthesisReference is
%      compared against the noise mask (mean + stdExcessFactor*std) and
%      against the max noise mask. If the reference PSD is entirely below 
%      the max noise mask, the window of the sourceSignal is considered 
%      noise, and the related filter for the window is null. In the other 
%      cases, the bins above the noise mask identify the spectrum to 
%      preserve. The filter mask is obtained by applying frequency and time 
%      clustering (frequency and time density masks) and smoothing in 
%      frequency and time by convolution with gaussian windows. The 
%      thresholds for the density masks are determined dynamically 
%      according to the level of the synthesisReference in respect of the 
%      noise level. 
% 
%   3- The smoothed filter is first decimated (see decimateFilter.m) and 
%      then antitransformed to obtain the time domain impulse response. The 
%      signal is filtered by forward linear convolution or forward/backward 
%      double linear convolution. Filtered windows are then overlapped, 
%      accounting for the leading and trailing tails in order to avoid 
%      artefacts.
%

% NOTES: 
% 
%   In the case of a simple filtering operation, the source signal 
%   (sourceSignal) is also the reference for the synthesis of the filter 
%   (synthesisReference). The latter does not need to be specified. In the
%   case of a signal-to-noise test for the filter performance, the 
%   synthesis reference can be loaded separately. Separating the source 
%   signal from the synthesis reference allows using the same filter to
%   process foreground and background signals separately. Note that the
%   reference signal and the noise mask used to process the foreground and 
%   background must be the same. See signalNoiseRatioTestUtility.m for an
%   application example using a synthetic signal. 
%



function [filteredSignal, sourceSignal, signalUnit, filterSummary] = noiseFilter(filterArgs)
    

    %LOAD ENVIRONMENT VARIABLES
    
    nameSpace
    spaceSettings


    %SETTINGS AND DEFAULT VALUES
    
    %signals
    if ~isfield(filterArgs, 'signals')
        
        error('Input signals not defined')
        
    else
        
        signals = filterArgs.signals;
    
        if isfield(signals, 'noiseMask') || ~isfield(signals, 'sourceSignal')
                    
        else    
            error('noise mask or input source signal not defined')
        end

        if ~isfield(signals, 'synthesisReference')
            signals.synthesisReference = [];
        end
        
    end

    pars.freqDensityLength_Hz         = [];
    pars.timeDensityLength_s          = [];
    pars.freqDensityThresholdMin      = [];
    pars.freqDensityThresholdMax      = [];
    pars.timeDensityThresholdMin      = [];
    pars.timeDensityThresholdMax      = [];
    pars.freqSmoothLength_Hz          = [];
    pars.timeSmoothLength_s           = [];
    pars.firFilterOrder               = [];
    pars.firFilterOrder_low           = [];
    pars.firFilterOrder_high          = [];
    pars.singleLinearConvolution      = [];
                
    if ~isfield(filterArgs, 'pars')
        filterArgs.pars = [];
    end
    if ~isfield(ss.filter, 'filter')
        ss.filter.filter = [];
    end

    pars = getSettings(pars, filterArgs.pars, ss.filter.filter, false);
   
    
    pictures.create = false;
    
    if isfield(filterArgs, 'pictures')
        
        if isfield(filterArgs.pictures, 'create')
            
            if filterArgs.pictures.create
            
                pictures.create = true;
                
                if isfield(filterArgs.pictures, 'frameTime_s')
                    pictures.frameTime_s = filterArgs.pictures.frameTime_s;
                else
                    pictures.frameTime_s = -1;
                end
                
                if isfield(filterArgs.pictures, 'outputFolder')
                    pictures.outputFolder = filterArgs.pictures.outputFolder;
                else
                    pictures.outputFolder = [];
                end     
                
                if isfield(filterArgs.pictures, 'format')
                    pictures.format = filterArgs.pictures.format;
                else
                    pictures.format = '';
                end     
                     
                if isfield(filterArgs.pictures, 'size')
                    pictures.size = filterArgs.pictures.size;
                else
                    pictures.size = '';
                end 
               
                if isfield(filterArgs.pictures, 'visible')
                    pictures.visible = filterArgs.pictures.visible;
                else
                    pictures.visible = true;
                end 
                    
                if isfield(filterArgs.pictures, 'plotWindow')
                    pictures.plotWindow = filterArgs.pictures.plotWindow;
                else
                    pictures.plotWindow = false;
                end 
            end
        end
    end
        
    
    summary.verboise = false;
    
    if isfield(filterArgs, 'summary')
        if isfield(filterArgs.summary, 'verboise')
            summary.verboise = filterArgs.summary.verboise;
        end
    end

            
    
    %LOAD SIGNALS
        
    %mask
    sampleRate      = signals.noiseMask.sampleRate;
    noiseMask       = signals.noiseMask.mask;
    maskLen         = length(signals.noiseMask.mask.mask);
    
    %window
    window          = signals.noiseMask.win;
    windowLen       = length(window.colaWindow);
    winStep         = (windowLen - window.overlap) / sampleRate;

    %Nyquist
    nf              = sampleRate/2;
    
    
    %source signal to filter
    try
        
        idx         = strfind(signals.sourceSignal, '\');
        idx         = idx(end);
        folder      = signals.sourceSignal(1:idx-1);
        fileName    = signals.sourceSignal(idx+1:end);
       
        [sourceSignal, signalUnit] = readAudio(folder, fileName, sampleRate);
        
    catch
        error('Signal audio file not found')
    end
    
    %reference
    if ~isempty(signals.synthesisReference)

        try

            idx                 = strfind(signals.synthesisReference, '\');
            idx                 = idx(end);
            folder              = signals.synthesisReference(1:idx-1);
            fileName            = signals.synthesisReference(idx+1:end);

            synthesisReference  = readAudio(folder, fileName, sampleRate);

        catch
            error('Reference audio file not found')
        end

        if length(sourceSignal) ~= length(synthesisReference)
            error('The source signal to filter must have the same length of the synthetis reference')
        end

    else
        synthesisReference   = [];
    end

    numberOfSamples = length(sourceSignal);
    
    if pictures.create == true && pictures.plotWindow == true
        
        plotWindowFunction(window.window, window.colaWindow, pictures);
    
    end
    
    
    
    %CHECK VALUES
    
    if  pars.freqDensityLength_Hz > 0 && pars.freqDensityLength_Hz <= nf
        
        densityFreqLen              = roundEven(pars.freqDensityLength_Hz/nf * (maskLen-1)) + 1;
        densityFreqLen              = max(1, densityFreqLen);
        pars.freqDensityLength_Hz   = (densityFreqLen - 1)/(maskLen-1) * nf;  
        
    else
        error('Density frequency range must be positive and shorter than Nyquist frequency')
    end

    if  pars.timeDensityLength_s > 0 && pars.timeDensityLength_s <= numberOfSamples/sampleRate
        
        densityTimeLen              = roundEven(pars.timeDensityLength_s/winStep) + 1;
        pars.timeDensityLength_s    = winStep * (densityTimeLen-1);
                
    else
        error('Density time range must be positive and shorter than the duration of the source signal')
    end

    if  pars.freqSmoothLength_Hz > 0 && pars.freqSmoothLength_Hz <= nf
        
        filterFreqSmoothLen         = round(pars.freqSmoothLength_Hz/nf * (maskLen-1)) + 1;
        filterFreqSmoothLen         = max(1, filterFreqSmoothLen);
        pars.freqSmoothLength_Hz    = (filterFreqSmoothLen - 1)/(maskLen-1) * nf;  
        
    else
        error('Filter smooth frequency range must be positive and shorter than Nyquist frequency')
    end

    if  pars.timeSmoothLength_s > 0 && pars.timeSmoothLength_s <= numberOfSamples/sampleRate
        
        filterTimeSmoothLen         = roundEven(pars.timeSmoothLength_s/winStep) + 1;
        pars.timeSmoothLength_s     = winStep * (filterTimeSmoothLen-1);
                
    else
        error('Filter smooth time range must be positive and shorter than the duration of the source signal')
    end
    
    
        
    %GET PSD
  
    [signalFft, ~, winCentralTime] =...
                stft(   sourceSignal,...
                        sampleRate,...
                        'Window',           window.colaWindow,...
                        'OverlapLength',    window.overlap,...
                        'FFTLength',        length(window.colaWindow));

    signalPsd = powerSpectralDensity(signalFft, sampleRate);
    
    if ~isempty(synthesisReference)

        [referenceFft, ~, ~] =...
                stft(   synthesisReference,...
                        sampleRate,...
                        'Window',           window.colaWindow,...
                        'OverlapLength',    window.overlap,...
                        'FFTLength',        length(window.colaWindow));     

        referencePsd = powerSpectralDensity(referenceFft, sampleRate);
    
    else
        
        referencePsd = signalPsd;
        
    end        
       
    numberOfWindows = size(signalPsd,2);
    
    
    
    %PROCESS WINDOWS
    
    %frame of the filter picture
    if pictures.create     
        
        if pictures.frameTime_s < 0 
            pictureFrameNumber = round(length(winCentralTime)/2);     
        elseif pictures.frameTime_s > winCentralTime(end)
            pictureFrameNumber = length(winCentralTime);
        else
            [~, frameIndex] = min(abs(pictures.frameTime_s - winCentralTime));
            pictureFrameNumber     = frameIndex;
        end
        
        pictures.frameTime_s = winCentralTime(pictureFrameNumber);
        
    end
    
        
    %prepare variables    
    densityFreqWin              = gausswin(densityFreqLen)';
    densityFreqWinSum           = sum(densityFreqWin);
    densityTimeWin              = gausswin(densityTimeLen)';
    densityTimeWinSum           = sum(densityTimeWin);
    
    densityFreqArrayLen         = maskLen + densityFreqLen-1;
    
   %densityFreqMatrix           = zeros(densityFreqArrayLen, densityFreqLen);
    densityTimeMatrix           = zeros(densityFreqArrayLen, densityTimeLen);
    
    timeSmoothPointer           = ceil(filterTimeSmoothLen/2);

    filtersSmoothTimeMatrix     = zeros(maskLen, filterTimeSmoothLen);
    densityThresholdCoeffBuffer = zeros(numberOfWindows,1);
    
    filteredSignal              = zeros(numberOfSamples, 1);
    

    for windowCounter = 1 : numberOfWindows + floor(densityTimeLen/2) + floor(filterTimeSmoothLen/2)

        if windowCounter <= numberOfWindows
        
            referenceWindow 	= referencePsd(:, windowCounter);
 
            %determine the dynamic density thresholds
            snMaxdiff = sort(referenceWindow - noiseMask.mask, 'descend');
            
            %0.1 every 10dB above the noise floor
            densityThresholdCoeff = snMaxdiff(round(maskLen*0.1))/10; 
            densityThresholdCoeff = max(0, densityThresholdCoeff);
            densityThresholdCoeff = min(1, densityThresholdCoeff);
            
            densityThresholdCoeffBuffer(windowCounter) = densityThresholdCoeff;
            
            densityFreqThreshold = pars.freqDensityThresholdMin +...
                                   densityThresholdCoeff*(pars.freqDensityThresholdMax - pars.freqDensityThresholdMin);
                                             
            %shift for a new item
            densityTimeMatrix 	= circshift(densityTimeMatrix, 1, 2);

            %boolean filters
            windowMaxFilterMask	= referenceWindow > noiseMask.max;   

            %compare masks
            signalExcessLevel 	= sum(windowMaxFilterMask)/length(noiseMask.max);

            if signalExcessLevel == 0

                densityTimeMatrix(:,1) = zeros(densityFreqArrayLen, 1);

            else

                %find the density of bins around each frequency bin
                densityFreqWindow = referenceWindow - noiseMask.mask;

                densityFreqWindow(densityFreqWindow < 0)    = 0;
                boolDensityFreqWindow                       = densityFreqWindow > 0;

                extendedBoolDensityWindow = [ zeros(floor(densityFreqLen/2), 1) + boolDensityFreqWindow(1);...  
                                              boolDensityFreqWindow;...
                                              zeros(floor(densityFreqLen/2), 1) + boolDensityFreqWindow(end)];

                boolDensityFreqMatrix      = repmat(extendedBoolDensityWindow, 1, densityFreqLen);              

                for shiftCounter = -floor(densityFreqLen/2) : floor(densityFreqLen/2)

                    boolDensityFreqMatrix(:, shiftCounter+ceil(densityFreqLen/2)) =...
                            circshift(boolDensityFreqMatrix(:, shiftCounter+ceil(densityFreqLen/2)), shiftCounter, 1);         

                end

                %weight the neighbours using a gaussian window
                densityFreqMatrix = boolDensityFreqMatrix .* repmat(densityFreqWin, densityFreqArrayLen, 1);   

                densityFreqArray  = sum(densityFreqMatrix, 2)/densityFreqWinSum > densityFreqThreshold;
         
                densityTimeMatrix(:,1) = densityFreqArray;

            end
            
        else
         
            densityTimeMatrix(:,1) = zeros(densityFreqArrayLen, 1);
            
        end
        
        if windowCounter > floor(densityTimeLen/2)
        
            %recall the density threshold coefficient
            if windowCounter - floor(densityTimeLen/2) <= numberOfWindows
                densityThresholdCoeff = densityThresholdCoeffBuffer(windowCounter - floor(densityTimeLen/2));
            else
                densityThresholdCoeff = densityThresholdCoeffBuffer(end);
            end
            
            densityTimeThreshold  = pars.timeDensityThresholdMin +...
                                    densityThresholdCoeff*(pars.timeDensityThresholdMin - pars.timeDensityThresholdMax);  
                  
            %weight the time neighbours using a gaussian window and sum
            timeDensityArray                = sum(densityTimeMatrix .* repmat(densityTimeWin, densityFreqArrayLen, 1), 2); 
            timeDensityArray                = timeDensityArray/densityTimeWinSum > densityTimeThreshold;
            %timeDensityArray                = timeDensityArray/densityTimeLen > densityTimeThreshold;

            %trim the ends
            densityFilterMask                = timeDensityArray(ceil(densityFreqLen/2) : end-floor(densityFreqLen/2));  
  
            %scale according to the intensity of the noise
            scaledFilterMask                = densityFilterMask .* (1 - (noiseMask.mask - mean(noiseMask.mask))/(max(noiseMask.mask) - mean(noiseMask.mask)));
            
            %frequency smooth filter
            freqSmoothedFilterMask          = smoothdata(scaledFilterMask,'gaussian',filterFreqSmoothLen);
            
            %nomalise
            maxSmoothedFilterMask           = max(max(freqSmoothedFilterMask),1);
            normalisedFilterMask            = freqSmoothedFilterMask/maxSmoothedFilterMask;

            %save filter mask   
            filtersSmoothTimeMatrix         = circshift(filtersSmoothTimeMatrix, 1, 2);
            filtersSmoothTimeMatrix(:,1)    = normalisedFilterMask;
            
           
            %time smooth filter 
            timeFreqSmoothedFilterMask      = smoothdata(filtersSmoothTimeMatrix, 2, 'gaussian', filterTimeSmoothLen);
            timeFreqSmoothedFilterMask      = timeFreqSmoothedFilterMask(:, timeSmoothPointer);


            %decimate filter mask
            [pars.firFilterOrder, decimatedFilterPulseResponse] = decimateFilter(pars.firFilterOrder,  timeFreqSmoothedFilterMask);

        end
            
        %account for the shift introduced by the time density and the time smooth
        windowShiftedPointer = windowCounter - floor(densityTimeLen/2) - floor(filterTimeSmoothLen/2);
        
        if windowShiftedPointer > 0
        
            %extract the time window
            startSample         = round(winCentralTime(windowShiftedPointer) * sampleRate - floor(windowLen/2) + 1);
            endSample           = startSample + windowLen - 1;
            timeWindow          = sourceSignal(startSample:endSample);
            timeWindow          = timeWindow .* window.colaWindow;

            overlapStartSample  = startSample;
            overlapEndSample    = endSample;

            %in sigle forward convolution the window extends only to the right 
            if pars.singleLinearConvolution == true

                filteredWindow = conv(decimatedFilterPulseResponse, timeWindow); 

                if windowShiftedPointer == numberOfWindows

                    filteredWindow = filteredWindow(1:windowLen);

                else

                    overlapEndSample = overlapEndSample + pars.firFilterOrder;

                end

            %in double forward/backward convolution the window extends both to the left and to the right
            else 

                linConv_I               = conv(decimatedFilterPulseResponse, timeWindow);
                revLinConv              = linConv_I(end:-1:1);
                linConv_II              = conv(decimatedFilterPulseResponse, revLinConv);
                filteredWindow          = linConv_II(end:-1:1);

                if windowShiftedPointer == 1

                    filteredWindow      = filteredWindow(pars.firFilterOrder+1:end);
                    overlapEndSample    = endSample + pars.firFilterOrder;

                elseif windowShiftedPointer == numberOfWindows

                    filteredWindow      = filteredWindow(1:windowLen + pars.firFilterOrder);
                    overlapStartSample  = startSample - pars.firFilterOrder;

                else

                    overlapStartSample  = startSample - pars.firFilterOrder;
                    overlapEndSample    = endSample + pars.firFilterOrder;

                end
            end

            %overlap and add
            filteredSignal(overlapStartSample:overlapEndSample) = filteredSignal(overlapStartSample:overlapEndSample) + filteredWindow;

            %create a picture of the filter for the frame required
            if pictures.create == true && pictureFrameNumber == windowShiftedPointer

                decimatedFilter         = fft(decimatedFilterPulseResponse);
                decimatedFilter         = decimatedFilter(1:floor((pars.firFilterOrder+1)/2)+1);
                decimatedFilter         = abs(decimatedFilter);

                [pars.firFilterOrder_low,  decimatedFilter_low] = decimateFilter(pars.firFilterOrder_low,  timeFreqSmoothedFilterMask);
                decimatedFilter_low    	= fft(decimatedFilter_low);
                decimatedFilter_low    	= decimatedFilter_low(1:ceil(pars.firFilterOrder_low/2)+1);
                decimatedFilter_low    	= abs(decimatedFilter_low);

                [pars.firFilterOrder_high, decimatedFilter_high] = decimateFilter(pars.firFilterOrder_high, timeFreqSmoothedFilterMask);
                decimatedFilter_high    = fft(decimatedFilter_high);
                decimatedFilter_high    = decimatedFilter_high(1:ceil(pars.firFilterOrder_high/2)+1);
                decimatedFilter_high    = abs(decimatedFilter_high);


                filters.desired.modulus         = timeFreqSmoothedFilterMask;
                filters.decimated.modulus       = decimatedFilter;
                filters.decimated.order         = pars.firFilterOrder;
                filters.decimated_low.modulus   = decimatedFilter_low;
                filters.decimated_low.order     = pars.firFilterOrder_low;
                filters.decimated_high.modulus  = decimatedFilter_high;
                filters.decimated_high.order    = pars.firFilterOrder_high;

                plotFilters(winCentralTime(windowShiftedPointer),...
                            noiseMask,...
                            sampleRate,...
                            referencePsd(:,windowShiftedPointer),...
                            filters,...
                            pictures);
            end
        
        end
        
    end 
    
    %output string
    
    filterSummary = createFilterSummary(signals, pars);
    
    if pictures.create || summary.verboise
    
        outputFilterSummary = filterSummary;
        outputFilterSummary = strcat(outputFilterSummary, sprintf("Figures output time: %.3f (s)\n", pictures.frameTime_s));
                
        %picture of the filtered signal and filter summary
        if pictures.create 
            
            plotSignals(sourceSignal, filteredSignal,  sampleRate, signalUnit,  pictures);

            if ~isempty(pictures.outputFolder)
                
                summaryFileId = fopen([pictures.outputFolder, '\', ns.filter.filterSummaryFileName, '.txt'], 'w');
                fprintf(summaryFileId,'%s', outputFilterSummary);
                fclose(summaryFileId);
                
            end
            
        end

        %prompt filter summary
        if summary.verboise
            
            fprintf(strrep(outputFilterSummary,'\','\\'));

        end
        
    end
        
end
        
        


%FUNCTIONS


% Calculate the single-sided PSD
%
% Inputs: 
%         fftArray: double-sided fast Fourier transform
%         sampleRate: sample rate
%
% Output: 
%         psdHalfDb: single-sided PSD

function psdHalfDb = powerSpectralDensity(fftArray, sampleRate)

    numberOfSamplesPerWindow = size(fftArray, 1);

    psdFull                 = 1/(sampleRate * numberOfSamplesPerWindow) * abs(fftArray).^2;
    psdHalf                 = psdFull(floor(numberOfSamplesPerWindow/2):end, :);    %take the second half including null and Nyquist frequency
    psdHalf(2:end-1)        = 2*psdHalf(2:end-1);                                   %add the second half without freq 0 e nyquist freq
    psdHalfDb               = 10*log10(psdHalf);

end



% Plot the filter model for the desired window
%
% Inputs: 
%         windowTime: central time of the window
%         noiseMask: noise mask
%         sampleRate: sample rate
%         referencePsd: PSD of the synthesis reference
%         filters: filter spectrums
%         pictures: picture properties

function plotFilters(windowTime, noiseMask, sampleRate, referencePsd, filters, pictures) 
    
    nameSpace

    settings.format         = pictures.format;
    settings.size           = pictures.size;
    plotSettings            = getPlotSettings(settings);
    
    
    nf                      = sampleRate/2; 
    frequencyBase           = linspace(0, nf/1000, length(noiseMask.mask));
    
    windowTimeString        = sprintf('Window center time:%.2f s', windowTime);
    filterModelFigure       = figure('name', windowTimeString, 'visible', pictures.visible);
    
    
    %mask subplot 
    noiseAveragePsdDb       = mean(noiseMask.mean);
    
    subplot(2,1,1)
    
    hold on 
    plot(frequencyBase, noiseMask.max     - noiseAveragePsdDb,...
            'LineWidth', plotSettings.lines.thin, 'Color', getColour('red'))
    plot(frequencyBase, noiseMask.mean    - noiseAveragePsdDb,...
            'LineWidth', plotSettings.lines.thin, 'Color', getColour('yellow'))
    plot(frequencyBase, noiseMask.mask    - noiseAveragePsdDb,...
            'LineWidth', plotSettings.lines.thin, 'Color', getColour('green'))
    plot(frequencyBase, referencePsd      - noiseAveragePsdDb,...
            'LineWidth', plotSettings.lines.thin, 'Color', getColour('blue'))
    hold off

    set(gca,'fontsize', plotSettings.labels.tickFontSize);
    
    yMin = min([(referencePsd - noiseAveragePsdDb)', (noiseMask.max - noiseAveragePsdDb)']);
    yMax = max([(referencePsd - noiseAveragePsdDb)', (noiseMask.max - noiseAveragePsdDb)']);  
    xlim([0, nf/1000]) 
    ylim([yMin, 1.3*yMax]) 
    
    xlabel('$Frequency~[kHz]$',...
            'FontSize', plotSettings.labels.fontSize, 'interpreter','latex')
    ylabel('$Amplitude~[dB]~~Re \overline{N_0}$',...
            'FontSize', plotSettings.labels.fontSize, 'interpreter','latex')

    lgdItemsM{1} = 'Max noise: N_{max}'; 
    lgdItemsM{2} = 'Mean noise: N_0';
    lgdItemsM{3} = sprintf('Mask: N_0 + %1.1f', noiseMask.stdExcessFactor);  lgdItemsM{3} = [lgdItemsM{3} '\sigma_N']; 

    [lgd,objh]              = legend(lgdItemsM, 'FontSize', plotSettings.legend.fontSize);    
    lgd.Location            = 'NorthEast';
    lgd.NumColumns          = 2;
    lgd.LineWidth           = plotSettings.legend.boxLineThickness;
    lgd.BoxFace.ColorType   = plotSettings.legend.boxColourType;
    lgd.BoxFace.ColorData   = plotSettings.legend.boxColourData;

    
    %filter subplot     
    filterFrequencyBase         = linspace(0, nf/1000, length(filters.decimated.modulus));
    filterFrequencyBase_low     = linspace(0, nf/1000, length(filters.decimated_low.modulus));
    filterFrequencyBase_high    = linspace(0, nf/1000, length(filters.decimated_high.modulus));
        
    subplot(2,1,2)
    
    hold on
    plot(frequencyBase,             filters.desired.modulus,...
                                    'LineWidth', plotSettings.lines.thin, 'Color', getColour('black'))
    plot(filterFrequencyBase,       filters.decimated.modulus,...       
            '-o' , 'MarkerSize',3,  'LineWidth', plotSettings.lines.thin, 'Color', getColour('green'))
    plot(filterFrequencyBase_low,   filters.decimated_low.modulus,....  
            '-s' , 'MarkerSize',3,  'LineWidth', plotSettings.lines.thin, 'Color', getColour('granate'))
    plot(filterFrequencyBase_high,  filters.decimated_high.modulus,...  
            '-x' , 'MarkerSize',3,  'LineWidth', plotSettings.lines.thin, 'Color', getColour('orange'))
    hold off

    set(gca,'fontsize', plotSettings.labels.tickFontSize);
    
    xlim([0, nf/1000])
    ylim([-0.05, 1.05]) 
    
    xlabel('$Frequency~[kHz]$',...
            'FontSize', plotSettings.labels.fontSize, 'interpreter','latex')
    ylabel('$|H(f)|$',...
            'FontSize', plotSettings.labels.fontSize, 'interpreter','latex')
        

    lgdItemsF{1} = 'Desired filter'; 
    lgdItemsF{2} = sprintf('FIR order %d', filters.decimated.order);
    lgdItemsF{3} = sprintf('FIR order %d', filters.decimated_low.order);
    lgdItemsF{4} = sprintf('FIR order %d', filters.decimated_high.order);

    [lgd,objh]              = legend(lgdItemsF, 'FontSize', plotSettings.legend.fontSize);    
    lgd.Location            = 'NorthEast';
    lgd.NumColumns          = 1;
    lgd.LineWidth           = plotSettings.legend.boxLineThickness;
    lgd.BoxFace.ColorType   = plotSettings.legend.boxColourType;
    lgd.BoxFace.ColorData   = plotSettings.legend.boxColourData;
        
    
    %set size
    set(filterModelFigure,...
        'Units',        'centimeters',...
        'Position',     plotSettings.position,...
        'PaperUnits',   'centimeters',...
        'PaperSize',    plotSettings.paperSize);

    
    %save pic
    if isempty(pictures.outputFolder)
        
        pause(5);
        
    else
        
        if ~exist(pictures.outputFolder, 'dir')
            mkdir(pictures.outputFolder);
        end

       
        fullFileName = [pictures.outputFolder, '\', ns.filter.filterModelFigureName, ' ', sprintf('(window center %.2f s)', windowTime), '.pdf'];   
        saveas(filterModelFigure, fullFileName);
        
    end
       
    close(filterModelFigure)
    
end



% Plot source and filtered signals
%
% Inputs: 
%         sourceSignal: source signal
%         filteredSignal: filtered signal
%         sampleRate: sample rate
%         signalUnit: signal physical unit (from source signal)
%         pictures: picture properties

function plotSignals(sourceSignal, filteredSignal,  sampleRate, signalUnit,  pictures)
    
    nameSpace

    settings.format         = pictures.format;
    settings.size           = pictures.size;
    plotSettings            = getPlotSettings(settings);
    
    
    numberOfSamples         = length(sourceSignal);
    duration                = length(sourceSignal)/sampleRate;
    t                       = linspace(0, duration, numberOfSamples);
    
    sourceVsFilteredFigure  = figure('name', 'source signal vs filtered signal', 'visible', pictures.visible);
    
 
    %source signal
    subplot(2,1,1)

    plot(t, sourceSignal, 'LineWidth', plotSettings.lines.thin)
        
    set(gca,'fontsize', plotSettings.labels.tickFontSize);
    
    xlim([0 t(end)])
    try
    ylim([-1.1 * max(abs(sourceSignal))  1.1 * max(abs(sourceSignal))])
    catch
    ylim([-1 1])   
    end
        
    xlabel('$Time~[s]$',...
                'fontsize', plotSettings.labels.fontSize, 'interpreter','latex');
    ylabel(['$Amplitude~', signalUnit, '$'],...
                'fontsize', plotSettings.labels.fontSize, 'interpreter','latex');
          
    lgdItemsM{1}            = 'Noisy'; 

    [lgd,objh]              = legend(lgdItemsM, 'FontSize', plotSettings.legend.fontSize);    
    lgd.Location            = 'NorthEast';
    lgd.NumColumns          = 1;
    lgd.LineWidth           = plotSettings.legend.boxLineThickness; 
    lgd.BoxFace.ColorType   = plotSettings.legend.boxColourType;
    lgd.BoxFace.ColorData   = plotSettings.legend.boxColourData;
            
            
    %filtered signal
    subplot(2,1,2)
    
    plot(t, filteredSignal, 'LineWidth', plotSettings.lines.thin)
    
    set(gca,'fontsize', plotSettings.labels.tickFontSize);
    
    xlim([0 t(end)]);
    try
    ylim([-1.1 * max(abs(filteredSignal))  1.1 * max(abs(filteredSignal))]);
    catch
    ylim([-1 1])  
    end
    
    xlabel('$Time~[s]$',...
                'fontsize', plotSettings.labels.fontSize, 'interpreter','latex');
    ylabel(['$Amplitude~', signalUnit, '$'],...
                'fontsize', plotSettings.labels.fontSize, 'interpreter','latex');

    lgdItemsM{1}            = 'Filtered'; 

    [lgd,objh]              = legend(lgdItemsM, 'FontSize', plotSettings.legend.fontSize);    
    lgd.Location            = 'NorthEast';
    lgd.NumColumns          = 1;
    lgd.LineWidth           = plotSettings.legend.boxLineThickness;
    lgd.BoxFace.ColorType   = plotSettings.legend.boxColourType;
    lgd.BoxFace.ColorData   = plotSettings.legend.boxColourData;
      
    
    %set size
    set(sourceVsFilteredFigure,...
        'Units',        'centimeters',...
        'Position',     plotSettings.position,...
        'PaperUnits',   'centimeters',...
        'PaperSize',    plotSettings.paperSize);
            
    %save
    if isempty(pictures.outputFolder)
        
        pause(5);
        
    else
        
        if ~exist(pictures.outputFolder, 'dir')
            mkdir(pictures.outputFolder);
        end

        fullFileName = [pictures.outputFolder, '\', ns.filter.sourceVsFilteredPicFileName, '.pdf'];   
        saveas(sourceVsFilteredFigure, fullFileName);
        
    end

    close(sourceVsFilteredFigure)
    
end
    
 

% Plot window function imported from getNoiseFrequencyMask.m
%
% Inputs: 
%         windowFunc: windowing function (Hann)
%         modifiedWindow: modified windowing for the COLA condition (depends on the overlap)
%         pictures: picture properties

function plotWindowFunction(windowFunc, modifiedWindow, pictures)
    
    nameSpace

    settings.format         = pictures.format;
    settings.size           = pictures.size;
    plotSettings            = getPlotSettings(settings);
    
    windowFunctionFigure   = figure('name', 'Window function', 'visible', pictures.visible);

    hold on
    plot(windowFunc,...
            '--', 'LineWidth', plotSettings.lines.thick, 'Color', getColour('black'))
    plot(modifiedWindow,...
            ':', 'LineWidth',  plotSettings.lines.thick, 'Color', getColour('red'))
    hold off
     
    xlim([0 length(windowFunc)])
    ylim([0  1.1])
    
    set(gca,'fontsize', plotSettings.labels.tickFontSize);
    xlabel('$Samples$',  'fontsize', plotSettings.labels.fontSize, 'interpreter','latex');
    ylabel('$Amplitude$','fontsize', plotSettings.labels.fontSize, 'interpreter','latex');
          
    lgdItems{1}            = 'NOLA'; 
    lgdItems{2}            = 'COLA'; 

    [lgd,objh]              = legend(lgdItems, 'FontSize', plotSettings.legend.fontSize);    
    lgd.Location            = 'south';
    lgd.NumColumns          = 1;
    lgd.LineWidth           = plotSettings.legend.boxLineThickness ;
    lgd.BoxFace.ColorType   = plotSettings.legend.boxColourType;
    lgd.BoxFace.ColorData   = plotSettings.legend.boxColourData;
    
    
    %set size
    set(windowFunctionFigure,...
        'Units',        'centimeters',...
        'Position',     plotSettings.position ,...
        'PaperUnits',   'centimeters',...
        'PaperSize',    plotSettings.paperSize);
            
    %save
    if isempty(pictures.outputFolder)
        
        pause(5);
        
    else
        
        if ~exist(pictures.outputFolder, 'dir')
            mkdir(pictures.outputFolder);
        end

        fullFileName = [pictures.outputFolder, '\', ns.filter.windowPicFileName, '.pdf'];   
        saveas(windowFunctionFigure, fullFileName);
        
    end

    close(windowFunctionFigure)
    
end

  

% Get filter summary string
%
% Inputs: 
%         signals: path of the signals
%         pars: filter parameters
%
% Output: 
%         filterSummary: filter summary string

function filterSummary = createFilterSummary(signals, pars)
    
    filterSummary = [];

    %signals
    if ~isempty(signals.synthesisReference)
        synthesisReferenceFile = signals.synthesisReference;
    else
        synthesisReferenceFile = signals.sourceSignal;
    end
            
    filterSummary = strcat( filterSummary,...
                            sprintf("AUDIO FILES:\n"),...
                            sprintf("Source signal:              %s\n",  signals.sourceSignal),...
                            sprintf("Synthesis reference signal: %s\n\n", synthesisReferenceFile));
                        
    %mask
    filterSummary = strcat( filterSummary,...
                            sprintf("MASK SETTINGS:\n"),...
                            signals.noiseMask.outputString,...
                            sprintf(" \n"));

    %filter
    stringFiterSettings          	= sprintf("FILTER SETTINGS:\n");
    sampleRateString           	    = sprintf("Sample rate: %d (samples per second)\n", signals.noiseMask.sampleRate);
    stringDensityFreq_Hz           	= sprintf("Frequency density interval length: %.1f (Hz)\n", pars.freqDensityLength_Hz);
    stringDensityTime_s            	= sprintf("Time density interval length: %.1f (s) \n", pars.timeDensityLength_s);
    stringDensityFreqThresholdMin  	= sprintf("Frequency density threshold min:  %.2f\n", pars.freqDensityThresholdMin);
    stringDensityFreqThresholdMax  	= sprintf("Frequency density threshold max:  %.2f\n", pars.freqDensityThresholdMax);
	stringDensityTimeThresholdMin 	= sprintf("Time density threshold min:  %.2f\n", pars.timeDensityThresholdMin);
	stringDensityTimeThresholdMax	= sprintf("Time density threshold max:  %.2f\n", pars.timeDensityThresholdMax);
	stringFilterFreqSmooth_Hz     	= sprintf("Frequency smooth interval length: %.1f (Hz)\n", pars.freqSmoothLength_Hz);
	stringFilterTimeSmooth_s     	= sprintf("Time smooth interval length: %.3f (s)\n", pars.timeSmoothLength_s);
	stringFirFilterOrder            = sprintf("Filter order: %d\n", pars.firFilterOrder);
	
    if pars.singleLinearConvolution
        stringSingleLinearConvolution     = sprintf("Convolution: Single forward\n");
    else
        stringSingleLinearConvolution     = sprintf("Convolution: Double forward/backward\n");  
    end
    
    filterSummary = strcat( filterSummary,...
                            stringFiterSettings,...  
                            sampleRateString,...
                            stringDensityFreq_Hz,...              
                            stringDensityTime_s,...               
                            stringDensityFreqThresholdMin,...      
                            stringDensityFreqThresholdMax,...      
                            stringDensityTimeThresholdMin,... 
                            stringDensityTimeThresholdMax,...   
                            stringFilterFreqSmooth_Hz,...         
                            stringFilterTimeSmooth_s,...    
                            stringFirFilterOrder,...
                            stringSingleLinearConvolution);
            
end



