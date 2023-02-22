% SIMPLE DEMO USING THE NOISE FILTER

% COPYRIGHT: 
%
%   MIT License
%   © 2023 Dario Chiantello <pipesoundlibrary@gmail.com>
%

clear
clc

%LOAD LIBRARY
loadLibray();
nameSpace
spaceSettings


%FOLDERS AND FILES

testFolder             = 'Test_Tap_1';
%testFolder              = 'Test_Disturbance_3';
noiseFileName           = 'background.wav';
inputFile              = [ns.repository.path, '\',  ns.soundbank.folderName, '\Taps\Tap_1\Event_Observations\Event_5\event_unfiltered.wav'];
%inputFile               = [ns.repository.path, '\',  ns.soundbank.folderName, '\Taps\Tap_1\Disturbance_Observations\Disturbance_3\disturbance_unfiltered.wav'];

noiseFolder             = [ns.repository.path, '\',  ns.soundbank.folderName, '\Taps\Tap_1\Background_Observations'];
outputFolder            = [ns.repository.path, '\',  'NoiseFilter', '\', testFolder];
inSignalResultFileName  = 'testFilter_NoisySignal.wav';   
outSignalResultFileName = 'testFilter_FilteredSignal.wav';  
resultSummaryFileName   = 'testFilter_Summary.txt';


%MASK AND FILTER
    maskArgs.inputFolder                                    = noiseFolder;
    maskArgs.inputFileName                                  = noiseFileName;
    %maskArgs.windowOverlap                                 = 35;   
    %maskArgs.disjointSet                                   = false;
filterArgs.signals.noiseMask                                = getNoiseFrequencyMask(maskArgs);
   filterArgs.signals.sourceSignal                          = inputFile;   
   filterArgs.pictures.create                               = true;
   filterArgs.pictures.plotWindow                           = true;
   filterArgs.pictures.outputFolder                         = outputFolder;
   filterArgs.pictures.format                               = '';
   filterArgs.summary.verboise                              = true;     
[filteredSignal, sourceSignal, signalUnit, filterSummary]   = noiseFilter(filterArgs);


%PLAY FILTERED SIGNAL
soundsc(filteredSignal, ss.defaultSampleRate);


%SAVE FILES
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

saveAudio(outputFolder, inSignalResultFileName,  sourceSignal,   filterArgs.signals.noiseMask.sampleRate, signalUnit, '');
saveAudio(outputFolder, outSignalResultFileName, filteredSignal, filterArgs.signals.noiseMask.sampleRate, signalUnit, '');


disp(newline())
disp('Done!')




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





