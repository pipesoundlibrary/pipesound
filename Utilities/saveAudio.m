% SAVE AUDIO FILE AND WITH THE DESIRED SAMPLE RATE

% COPYRIGHT: 
%
%   MIT License
%   © 2023 Dario Chiantello <pipesoundlibrary@gmail.com>
%

% INPUTS: 
%   
%   outputFolder: path of the output folder
%   signalFileName: name of the saved .wav or .mp4 audio file (with extension)
%   signal: audio signal
%   sampleRate: audi sample rate
%   unit: string indicating the physical unit for the metadata. E.g. ([PA])
%   fileId: stored in the title field of the file metadata
%

function saveAudio(outputFolder, signalFileName, signal, sampleRate, unit, fileId)
    
    try
        
        dotPos          = strfind(signalFileName, '.');
        dotPos          = dotPos(end)+1;
        fileExtension   = signalFileName(dotPos:end);
        
    catch
        
        error('Invalid file name')
        
    end
        
    if strcmpi(fileExtension, 'wav') && strcmpi(fileExtension, 'mp4') 
    
       error('Save .wav or .mp4 files')
        
    end

    scaleFactor   	= max(abs(signal));
    normSignal    	= signal/scaleFactor;

    fileTitle       = fileId;
    fileComment     = setUnitAndScale(scaleFactor, unit);
    
    outputFile      = [outputFolder, '\', signalFileName];
    audiowrite(outputFile, normSignal, sampleRate, 'Comment', fileComment, 'Title', fileTitle);
    
end


