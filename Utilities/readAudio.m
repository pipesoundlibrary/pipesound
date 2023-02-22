% READ AUDIO FILE AND SET THE DESIRED SAMPLE RATE

% COPYRIGHT: 
%
%   MIT License
%   © 2023 Dario Chiantello <pipesoundlibrary@gmail.com>
%

% INPUTS: 
%   
%   signalFolder: audio signal folder
%   signalFileName: signal file name in .wav or .mp4 format (with extension)
%   sampleRate: resample rate of the audio file
% 

% OUTPUT: 
%   
%   signal: imported signal resampled @sampleRate
%   unit: physical unit string obtained from the file metadata (e.g. [PA])
% 

function [signal, unit] = readAudio(signalFolder, signalFileName, sampleRate)

    dotPos          = strfind(signalFileName, '.');
    dotPos          = dotPos(end)+1;
    fileExtension   = signalFileName(dotPos:end);

    if strcmpi(fileExtension, 'wav') && strcmpi(fileExtension, 'mp4') 
    
       error('Select .wav or .mp4 files')
        
    end
        
    signalFile = [signalFolder, '\', signalFileName];
    
    [audio, audioSampleRate]   = audioread(signalFile);

    audioMetadata              = audioinfo(signalFile);
    [audioScaleFactor, unit]   = getUnitAndScale(audioMetadata);
    
    if sampleRate ~= audioSampleRate
        [p,q]                 = rat(sampleRate / audioSampleRate);
        audio                 = resample(audio,p,q);
    end
    
    signal = audioScaleFactor * audio; 
    
    
    