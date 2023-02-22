% MODE EXTRACTOR FOR THE RESULTS OBTAINED FROM KZSOLVER.

% COPYRIGHT: 
%
%   MIT License
%   © 2023 Dario Chiantello <pipesoundlibrary@gmail.com>
%

% DESCRIPTION:
% 
%   The extractor processes the results obtained from kzsolver.m to 
%   separate and sort the numerical values of the modes stored in  
%   solverResults_XXy.mat
% 

% INPUTS: 
%   
%   Input data are obtained from the intermediate processing result file 
%   solverResults_XXy.mat generated by the solver. Input settings can be 
%   provided through an external setting file or obtained from 
%   spaceSettings.m as default values. The customised setting definition 
%   requires the file extractorSetup.m stored under 
%   [simulation folder]\SimulationSetup. Examples can be found in 
%   Examples\AcousticModels. The setting structs are reported below, along 
%   with the related description. Fields not specified are obtained from 
%   spaceSettings.m as default values. 
%
%
%   extractorSetup.m
%
%   -exs.settings.
%      nMin                          Min angular mode order n
%      nMax                          max angular mode order n (N)
%      evanProp                      0-only propagative, 1-only evanescent, 2-both
%      maxContinuityRadius_rad_m     Max distance radius in rad/m for consecutive points of the same branch.
%      differentialAverageSteps      Number of frequency steps for the average of the differential. Used to sort the branches and close the gaps.
%      minBranchLen                  Minimum length in frequency steps of a branch chunck.
%      initialGapMinK_rad_m          Max value in rad/m above which the initial gap can be closed.
%      initialGapMinFreq_Hz          Max value in Hz to close the initial gap.
% 
%   -exs.pictures.
%      fLim_Hz                       Plot frequency limits in Hz (e.g. [0, 50000]). Leave empty to use the full solver range.
%      kLim_rad_m                    Plot kz limits in rad/m (e.g. [0, 250]).       Leave empty to use the full solver range. 
%      colours                       Colors for the plotted modes. Leave empty for default settings.
%      format                        Picture format, see Utilities/getPlotSettings 
%      size                          Picture size, see Utilities/getPlotSettings  
%

% OUTPUTS: 
%   
%   - modes.mat: table containing the extracted modes. The table is
%     partially overwritten when a new set of modes is extracted. Placed in 
%     [simulation folder]\Modes.
%   - modes_XXy.pdf: representation of the modes in the frequency-kz plane.
%     Placed in [simulation folder]\Modes\Mode_XX\modeType
%   - cutFrequencies_XXy.xls: cut-on/off frequency table. Placed in 
%     [simulation folder]\Modes\Mode_XX\modeType
% 

% PROCESSING STEPS: 
% 
%   1- The solution obtained through the solver is sorted by frequency 
%      step.
% 
%   2- Points belonging to the same branch and the same mode are sorted by 
%      distance using a differential prediction.
% 
%   3- Gaps are closed in the solution by linear interpolation.
% 
%   4- Names of the radial mode orders are assigned manually. The mode data 
%      are collected and stored in modes.mat. Only the modes extracted are 
%      overwritten in the table. A cut-on/off excel table is also 
%      generated.
% 

% NOTES: 
% 
%   - No settings to edit in this file. Just run and select the simulation.
%   - To use customised settings, add an extractor setup file in the 
%     simulation folder. Examples of extractorSetup.m can be found in 
%     Examples\AcousticModels
% 


clear
clc

%LOAD LIBRARY
loadLibray();
nameSpace
spaceSettings


%SELECT A SIMULATION FOLDER
simulationFoldersPath   = [ns.repository.path, '\', ns.acMod.folderName];
simulationFolders       = getFolderNames(simulationFoldersPath);

fprintf("Select a simulation folder:\n\n");

for simulationCounter = 1:length(simulationFolders)

    fprintf("%d > %s\n", simulationCounter-1, simulationFolders{simulationCounter});
    
end

prompt        = sprintf("\n>");
simulationId  = input(prompt, "s");

try
    
    simulationId = round(str2double(simulationId));
    
    if simulationId <0 || simulationId >= length(simulationFolders) || isnan(simulationId)
        simulationId = -1;   
    end
        
catch
    
    simulationId = -1;
    
end

if simulationId < 0
    error("Select a valid simulation")
else
    simulationFolderName = simulationFolders{simulationId+1};
    fprintf("Simulation selected: %s\n",simulationFolderName );
end





%LOAD SIMULATION FOLDER  
simulationFolderPath = [ns.repository.path, '\', ns.acMod.folderName,'\',simulationFolderName];
modeRootDirectory    = [simulationFolderPath,'\', ns.acMod.modes.folderName];

if ~exist(modeRootDirectory, 'dir')
   error('Run kzSolver before modeExtractor') 
end


%OVERWRITE
overwrite      = false;
prompt         = sprintf("Do you want to save and overwrite extractor result data? [Y/N] >\n");
overwriteAns   = input(prompt, "s");
if strcmp(overwriteAns, 'Y') ||  strcmp(overwriteAns, 'y')
    overwrite = true;
else
    warning("New data won't not be saved");
end


%LOAD EXTRACTOR SETTINGS  
simulationSetupFolderPath               = [simulationFolderPath,'\',ns.acMod.setup.folderName];
addpath(simulationSetupFolderPath);

if exist([simulationSetupFolderPath, '\', ns.acMod.setup.extractor.fileName, '.m'], 'file')
    
    exsIn = extractorSetup();
    
else
    
    wString = sprintf('Extractor setup file missing in: %s\n', simulationSetupFolderPath);
    wString = [wString, 'The extractor will use the default values defined in spaceSettins.m'];
    warning(wString)    
    exsIn = [];
    
end

%extractor settings
exs.settings.nMin                       = []; 
exs.settings.nMax                      	= [];                                      
exs.settings.evanProp                 	= [];                                      
exs.settings.maxContinuityRadius_rad_m 	= [];                                       
exs.settings.differentialAverageSteps 	= []; 
exs.settings.minBranchLen              	= [];                                        
exs.settings.initialGapMinK_rad_m     	= [];                                    
exs.settings.initialGapMinFreq_Hz     	= [];                            

if ~isfield(exsIn, 'settings')
    exsIn.settings = [];
end
if ~isfield(ss.acMod, 'extractor')    
    ss.acMod.extractor = [];
end
exs.settings = getSettings(exs.settings, exsIn.settings, ss.acMod.extractor, false);

%figure settings
exs.pictures.fLim_Hz                    = [];             
exs.pictures.kLim_rad_m                 = [];             
exs.pictures.colours                    = {};    
exs.pictures.format                     = '';      
exs.pictures.size                       = '';     

if isfield(exsIn, 'pictures')

    if isfield(exsIn.pictures, 'fLim_Hz')
        exs.pictures.fLim_Hz = exsIn.pictures.fLim_Hz; 
    end
    
    if isfield(exsIn.pictures, 'kLim_rad_m')
        exs.pictures.kLim_rad_m = exsIn.pictures.kLim_rad_m; 
    end
    
    if isfield(exsIn.pictures, 'colours')
        exs.pictures.colours = exsIn.pictures.colours; 
    end
    
    if isfield(exsIn.pictures, 'format')
        exs.pictures.format = exsIn.pictures.format;
    end

    if isfield(exsIn.pictures, 'size')
        exs.pictures.size = exsIn.pictures.size;
    end
end



%output table n modes
nModesTable     = table();

%propagative/evanescent selection
if exs.settings.evanProp <= 1
    evanProp = exs.settings.evanProp;
else
    evanProp = 0:1;
end    

%loop over N group of modes
for n = exs.settings.nMin:exs.settings.nMax
    
    %check if the mode result folder exists
    modeFolderName  = [ns.acMod.modes.modeFolderRootName, '_', num2str(n), 'X'];
    modeFolder      = [modeRootDirectory, '\', modeFolderName];
    if ~exist(modeFolder, 'dir')
        warning('Result folder of mode %s missing.', modeFolderName) 
        continue
    end
    
    %prepare output tables
    nTable      = table(cell(1), 'VariableNames', {['n=',num2str(n)]});
    nModesTable = [nModesTable, nTable];
    
    modes       = table(cell(1),cell(1), 'VariableNames', {ns.acMod.modes.propagativeFolderName, ns.acMod.modes.evanescentFolderName});   
    nModesTable.(num2str(['n=',num2str(n)])) = modes;
    
  
    %evanescent/propagative
    for evanescent = evanProp
    
        %load solver results  
        if evanescent
            
            resultsFileName     = [ns.acMod.kzSolver.resultRootName, '_', num2str(n), 'Xe'];
            resultFolder        = [modeFolder, '\', ns.acMod.modes.evanescentFolderName];
            modeTypeString      = ns.acMod.modes.evanescentFolderName;

        else
            
            resultsFileName     = [ns.acMod.kzSolver.resultRootName, '_', num2str(n), 'Xp'];
            resultFolder        = [modeFolder, '\', ns.acMod.modes.propagativeFolderName];
            modeTypeString      = ns.acMod.modes.propagativeFolderName;

        end
        
        fprintf('\n**********\n')
        fprintf('%s mode n=%d\n\n', modeTypeString, n)
        
        
        resultsFile = [resultFolder '\', resultsFileName, '.mat' ];
        if ~exist(resultsFile, 'file')
            warning('Result file %s missing.', resultsFile) 
            continue
        end
        load(resultsFile);


                
        %organise points by frequency step
        separatedModes = zeros(1, length(fVector));
        
        for pointCounter = 1:size(compactModeMatrix,1)
        
            fStep       = find(compactModeMatrix(pointCounter,1) == fVector);
            nextEmpty   = sum(separatedModes(:,fStep) > 0) + 1;
            
            if nextEmpty > size(separatedModes,1)
            
                separatedModes = [separatedModes; zeros(1, length(fVector))];
                
            end
            
            separatedModes(nextEmpty, fStep) = compactModeMatrix(pointCounter, 2); 
            
        end
        
        
        
        %sort by distance using a differential prediction
        numberOfBranches = find(separatedModes(:,1)>0, 1, 'last');
        if isempty(numberOfBranches)
            numberOfBranches = 0;
        end
        
        differentialWeights = 0 : 1/exs.settings.differentialAverageSteps : 1;
        
        for fStep = 2:length(fVector)
        
            %current and previous values
            prevStepPoints      = separatedModes(:, fStep-1);
            currStepPoints      = separatedModes(:, fStep);
            
            %differentials
            if fStep > exs.settings.differentialAverageSteps+1
                
                prevRange = separatedModes(:, fStep-1-exs.settings.differentialAverageSteps:fStep-1);
                                
                for lineCounter = 1:size(prevRange,1)
                
                    selectedLine = prevRange(lineCounter, :);
                    
                    samples      = find(selectedLine>0);
                    sampleVals   = selectedLine(samples);
                    
                    if length(samples) >=2
                        prevRange(lineCounter, :) = interp1(samples, sampleVals, 1:exs.settings.differentialAverageSteps+1,'linear','extrap');
                    else
                        valuePos = prevRange(lineCounter, :)>0;
                        
                        if sum(valuePos > 0)  
                            prevRange(lineCounter, 1:end) = prevRange(lineCounter, valuePos);
                        end
                    end
                end
                
                differentials           = prevRange - circshift(prevRange, 1, 2);
                weightsMat              = repmat(differentialWeights,size(prevRange,1),1);               
                averagedDifferentials   = sum(weightsMat.*differentials,2)/sum(differentialWeights);
                
                predictedValues         = prevRange(:,end) + averagedDifferentials;
                                
            else
                
                predictedValues         = prevStepPoints;
                
            end
            
            %find the next by min distance
            newStepPoints = zeros(size(currStepPoints));
            
            ignoreZeroCurrStepPoints = currStepPoints;
            ignoreZeroCurrStepPoints(currStepPoints == 0)   = -10 * kVector(end); %just to ignore null values when calculating the min
            ignoreZeroPredictedValues = predictedValues;
            ignoreZeroPredictedValues(predictedValues == 0) = 10 * kVector(end);  %just to ignore null values when calculating the min
            
            [pointIds, minDiff] = findClosestNext(ignoreZeroPredictedValues, ignoreZeroCurrStepPoints, sum(currStepPoints>0));
            
            for pointJoinCounter = 1:size(pointIds,1)
                
                if minDiff(pointJoinCounter) < exs.settings.maxContinuityRadius_rad_m

                    newStepPoints(pointIds(pointJoinCounter,1))  = currStepPoints(pointIds(pointJoinCounter,2));
                    currStepPoints(pointIds(pointJoinCounter,2)) = 0;
                    
                end
            
            end
            
            
            %add the rest at the bottom           
            freeLines = length(newStepPoints) - numberOfBranches;
            newLines  = sum(currStepPoints>0) - freeLines;
            
            if newLines > 0
                
                newStepPoints  = [newStepPoints; zeros(newLines,1)];
                separatedModes = [separatedModes; zeros(newLines, length(separatedModes))];
                
            end
             
            currStepPointsIds = find(currStepPoints); 
            newStepPoints(numberOfBranches+1 : numberOfBranches + length(currStepPointsIds)) = currStepPoints(currStepPointsIds);

            
            
            %update
            separatedModes(:, fStep) = newStepPoints;
            
            if numberOfBranches < (find(newStepPoints>0,1,'last'))

                numberOfBranches = (find(newStepPoints>0,1,'last'));

            end 
            
        end
        
                
        %clean array
        separatedModes = separatedModes(sum(separatedModes>0, 2) > exs.settings.minBranchLen, :); 
        
        %close gaps 
        M = size(separatedModes,1);
        differentials = circshift(separatedModes,-1,2)- separatedModes;
        
        for m = 1:M
        
            valuePos                = separatedModes(m, :) > 0;
            valuePosLen             = find(valuePos, 1, 'last');
            
            firstPos                = find(valuePos, 1, 'first');
            valueFirstPos           = separatedModes(m, firstPos);
            
            [~, initialGapMinStep]  = min(abs(fVector - exs.settings.initialGapMinFreq_Hz));
             
            
            if firstPos > 1 && firstPos < initialGapMinStep && valueFirstPos > exs.settings.initialGapMinK_rad_m && valueFirstPos < 0.99*kVector(end) 
                gapStart = 1;
            else
                gapStart = 0;
            end
            
            gapEnd = 0;
            
            for posCounter = 2:valuePosLen
                       
                if valuePos(posCounter-1) == 1 && valuePos(posCounter) == 0
                    gapStart = posCounter;
                end

                if valuePos(posCounter-1) == 0 && valuePos(posCounter) == 1 && gapStart > 0
                    gapEnd = posCounter;
                end
                
                
                %fill gap
                if gapStart > 0 && gapEnd > 0
                
                    valueGapEnd     = separatedModes(m, gapEnd);
                    
                    if gapStart == 1                       
                        valueGapStart = valueGapEnd;
                    else
                        valueGapStart = separatedModes(m, gapStart-1);
                    end
                    
                    gapValues = interp1([gapStart-1, gapEnd], [valueGapStart, valueGapEnd], [gapStart-1 : gapEnd]);
                    
                    separatedModes(m, gapStart: gapEnd-1) = gapValues(2:end-1);
                        
                    gapStart = 0;
                    gapEnd   = 0;
                
                end
            end      
        end
        
        
        %assign mode names
        nodeNames       = cell(M,1);
        
        figureName      = sprintf('%s modes (n=%d). Names definition.', modeTypeString, n);
        modeNameFigure  = figure('name', figureName);
        
        for m = 1:M
            
            hold on
            for selection = 1:size(separatedModes,1)

                if selection == m
                    lineColour   = 'red'; 
                else
                     lineColour = 'black';
                end
                
                firstPoint = find(separatedModes(selection,:)>0, 1, 'first');
                lastPoint  = find(separatedModes(selection,:)>0, 1, 'last');
                
                plot(fVector(firstPoint:lastPoint), separatedModes(selection,firstPoint:lastPoint), 'Color', lineColour, 'LineWidth', 1.5);

            end
            hold off
            
            if evanescent
                defaultName = char(64 + m);
            else
                defaultName = num2str(m);
            end
            
            prompt      = sprintf("Please enter the name of the red mode. Press enter for the default name: %s) >\n", defaultName);
            mModeName   = input(prompt, "s");
            
            if isempty(mModeName)
                mModeName = defaultName;
            end
            
            nodeNames{m} = mModeName;
            
        end
        
        close (modeNameFigure)   
        
        
        
        %create mode table
        mModes = table();
        
        for m = 1:M
       
            frequency   = fVector';
            kz          = separatedModes(m,:)';
        
            frequency   = frequency(kz>0);
            kz          = kz(kz>0);
            
            values      = table(frequency, kz, 'VariableNames', {'frequency (Hz)', 'kz (rad/m)'});
            valTable    = table({values}, 'VariableNames', {['m=',nodeNames{m}]});
            
            mModes      = [mModes, valTable];
            
        end
        
        [~,sortedNamesIds]  = sort(mModes.Properties.VariableNames);
        mModes              = mModes(:,sortedNamesIds);
        
        prompt = sprintf('How many modes do you want to include? (Default: %d) >\n', M);
        try
            numberOfmodesSelected = round(str2double(input(prompt, "s")));
        catch 
            numberOfmodesSelected = M; 
        end
        if isnan(numberOfmodesSelected)
            numberOfmodesSelected = M; 
        end
        numberOfmodesSelected = min(M,numberOfmodesSelected);
        
        fprintf('Number of modes included: %d\n\n', numberOfmodesSelected);
        
        
        
        %add zero-cut from vector zeroCutFrequencies and create cut-on/off table
        fistValues  = zeros(M,3); %freq, kz, flag: -1 scan limits or same branch, +1 zero cut, 0 the others
        lastValues  = zeros(M,3); %freq, kz, flag: -1 scan limits or same branch, +1 zero cut, 0 the others
        cutOnFreq   = zeros(M,1); 
        cutOffFreq  = zeros(M,1); 
        
        
        %scan limits
        for m = 1:M
            
            fistValues(m,1:2) = mModes.(m){1}{1,:};
            lastValues(m,1:2) = mModes.(m){1}{end,:};
            
            if fistValues(m,1) == fVector(1) || fistValues(m,2) > 0.98 * kVector(end)
            
                fistValues(m,3) = -1;           %first value ok
                cutOnFreq(m)    = fVector(1); 	%cut-on at the first freq
                 
            end
                
            if lastValues(m,1) > 0.99 * fVector(end) || lastValues(m,2) > 0.99 * kVector(end)
            
                lastValues(m,3) = -1;       	%last value ok
                cutOffFreq(m)   = fVector(end); %cut-ff at the last freq
                
            end
        end
        
        
        %check if two modes belong to the same branch of the solution
        fProximityRadius = fVector(3) - fVector(1);                         % 2 frequency steps
        kProximityRadius = 5 * exs.settings.maxContinuityRadius_rad_m;      % k proximity (expand radius, bending point of the branch)
        
        for m1 = 1:M
            for m2 = m1+1:M
            
                if abs(fistValues(m1,1) - fistValues(m2,1)) <= fProximityRadius &&...   
                   abs(fistValues(m1,2) - fistValues(m2,2)) <= kProximityRadius     
                
                   fistValues(m1,3) = -1;
                   fistValues(m2,3) = -1;
                   cutOnFreq(m1)    = fistValues(m1,1);
                   cutOnFreq(m2)    = fistValues(m2,1);
                   
                end
                   
                if abs(lastValues(m1,1) - lastValues(m2,1)) <= fProximityRadius &&... 
                   abs(lastValues(m1,2) - lastValues(m2,2)) <= kProximityRadius    
                
                   lastValues(m1,3) = -1;
                   lastValues(m2,3) = -1;
                   cutOffFreq(m1)   = lastValues(m1,1);
                   cutOffFreq(m2)   = lastValues(m2,1);
                   
                end

            end
        end
        
        
        %add zero-cut and create cut-on/off tables
        fvFreq = fistValues(:,1);
        checkFirstValues = fistValues(:,3) == 0; 
        fvFreq(~checkFirstValues) = -fVector(end);
       
        lvFreq = lastValues(:,1);
        checkLastValues = lastValues(:,3) == 0; 
        lvFreq(~checkLastValues) = -fVector(end);
        
        numOfMins         = min(sum(checkFirstValues) + sum(checkLastValues), length(zeroCutFrequencies));
        [minIds, minDist] = findClosestNext([fvFreq; lvFreq], zeroCutFrequencies', numOfMins);
        
        for newCutPointCounter = 1:numOfMins
            
            if minIds(newCutPointCounter, 1) <= M

                cutOnFreq(minIds(newCutPointCounter, 1))         = zeroCutFrequencies(minIds(newCutPointCounter, 2));
                fistValues(minIds(newCutPointCounter, 1), 3)	 = 1;

            else

                cutOffFreq(minIds(newCutPointCounter, 1) - M)    = zeroCutFrequencies(minIds(newCutPointCounter, 2));
                lastValues(minIds(newCutPointCounter, 1) - M, 3) = 1;

            end
        end
        
        
        %add the cut-on/off values to the output mode table where necessary
        for m = 1:M
        
            if fistValues(m,3) > 0
            
                temp            = mModes.(m){1};
                newLine         = table(cutOnFreq(m), 0, 'VariableNames', temp.Properties.VariableNames);
                temp            = [newLine; temp];
                mModes.(m)(1)   = {temp};
             
            elseif fistValues(m,3) == 0
                
                cutOnFreq(m)    = fistValues(m,1);
                
            end
                
            if lastValues(m,3) > 0
                
                temp           	= mModes.(m){1};
                newLine       	= table(cutOffFreq(m), 0, 'VariableNames', temp.Properties.VariableNames);
                temp          	= [temp; newLine];
                mModes.(m)(1) 	= {temp};
                 
            elseif lastValues(m,3) == 0
                
                cutOffFreq(m)   = lastValues(m,1);
                
            end
        end
        
        
        
        %OUTPUT  
        
        mModesTable     = mModes(:,1:numberOfmodesSelected);
        cutOnVector     = cutOnFreq(1:numberOfmodesSelected,1);
        cutOffVector    = cutOffFreq(1:numberOfmodesSelected,1);
           
        %save cut-on/cut-off frequency table
        saveCutFrequencies(resultFolder, mModesTable, n, evanescent, cutOnVector, cutOffVector, overwrite);
        
        %print and save mode figures        
        printModeFigure(resultFolder, mModesTable, n, evanescent, fVector, kVector, exs.pictures, overwrite);
        
        %save data
        nModesTable.(['n=',num2str(n)]).(modeTypeString) = mModesTable;
        saveData(modeRootDirectory, nModesTable,  n, evanescent, overwrite);

        
    end %for evanescent/propagative
    
end %for n

fprintf('Done!\n')




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




% Calculates the mutually exclusive minimums for two arrays
%
% Inputs:  
%           currentVector: first vector of values
%           nextVector: second vector of values
%           numberOfNextPoints: number of minimum point to find
%
% Outputs: 
%           minIds: indexes of the closest items for both input vectors.
%                   One line stores the indexes of two coupled items, from 
%                   the first line for the closest to the last for the 
%                   furthest. 
%           minDist: distances ordered from the closest to the furthest  
%

function [minIds, minDist] = findClosestNext(currentVector, nextVector, numberOfNextPoints)
        
    currentLen          = length(currentVector);
    nextLen             = length(nextVector);

    numberOfNextPoints  = min([numberOfNextPoints currentLen nextLen]);
    
    currentValuesMat    = repmat(currentVector, 1, nextLen);
    nextValuesMat       = repmat(nextVector', currentLen, 1);
    diffMat             = abs(currentValuesMat - nextValuesMat);
    maxDiffMat          = max(max(diffMat))+1;
    
    minIds              = zeros(numberOfNextPoints, 2);
    minDist             = zeros(numberOfNextPoints, 1);
    
    for minCounter = 1:numberOfNextPoints

        [minVals, minRows]      = min(diffMat);
        [minVal,  minCol]       = min(minVals);
        minRow                  = minRows(minCol);

        minIds(minCounter, :)   = [minRow, minCol];
        minDist(minCounter)     = diffMat(minRow, minCol);
        
        diffMat(:, minCol)      = maxDiffMat;
        diffMat(minRow, :)      = maxDiffMat;
    end  
    
end




% Print a table containing the cut-on/off frequencies and save .xls
%
% Inputs:  
%           resultFolder: output folder
%           mModesTable: mode result table
%           n: mode order
%           evanescent: mode type (true/false)
%           cutOnFreq: cut-on frequency vector
%           cutOffFreq: cut-off frequency vector
%           overwrite: save and overwrite data (true/false)
%

function saveCutFrequencies(resultFolder, mModesTable, n, evanescent, cutOnFreq, cutOffFreq, overwrite)

    modeNameStrings = mModesTable.Properties.VariableNames;
    matrixFields    = [{''}; {'Cut-on (kHz)'}; {'Cut-off (kHz)'}];
    cutFrequencies  = [modeNameStrings; num2cell(cutOnFreq/1000)'; num2cell(cutOffFreq/1000)'];
    cutFrequencies  = [matrixFields, cutFrequencies];

    %print output
    if evanescent
        fprintf('Cut frequencies for evenescent modes n=%d\n', n);
    else
        fprintf('Cut frequencies for propagative modes n=%d\n', n);
    end
    cutFrequencies
    
    if overwrite
        
        nameSpace
        
        if evanescent
            modeId =  [num2str(n), 'Xe',];
        else
            modeId =  [num2str(n), 'Xp'];
        end
        
        fileName        = [resultFolder, '\', ns.acMod.extractor.cutFreqFileName, '_', modeId]; 
        xlswrite(fileName, cutFrequencies, 'Output','A1');
    
    end
    
end




% Print mode figure and save .pdf
%
% Inputs:  
%           outputDirectory: output folder
%           mModesTable: mode result table
%           n: mode order
%           evanescent: mode type (true/false)
%           fVector: frequency vector
%           kVector: kz vector
%           settings: figure settings (see extractor setup)
%           overwrite: save and overwrite data (true/false)
%

function printModeFigure(outputDirectory, mModesTable, n, evanescent, fVector, kVector, pictureSettings, overwrite)

    nameSpace

    plotSettings = getPlotSettings(pictureSettings);
            
    M = size(mModesTable,2);
    
    if isempty(pictureSettings.fLim_Hz)
        fLim = [fVector(1), fVector(end)];
    else
        fLim = pictureSettings.fLim;
    end

    if isempty(pictureSettings.kLim_rad_m)
        kLim = [kVector(1), kVector(end)];
    else
        kLim = pictureSettings.kLim;
    end

    if evanescent
        figureNameString = sprintf('Evenescent modes n=%d\n', n);
    else
        figureNameString = sprintf('Propagative modes n=%d\n', n);
    end
    
    modesFigure = figure('Name', figureNameString);
    
    lgd = mModesTable.Properties.VariableNames;   
    lgd = cellfun(@(x) sprintf('(%d,%s)', n, x(3:end)), lgd, 'UniformOutput', false);
    
    colours = [];
    if ~isempty(pictureSettings.colours)
    
        if evanescent 
            if isfield(pictureSettings.colours, 'evanescent')
                colours = pictureSettings.colours.evanescent;
            end
        else
            if isfield(pictureSettings.colours, 'propagative')
                colours = pictureSettings.colours.propagative;
            end
        end
        
    end
    
    hold on
    for m = 1:M

        values = mModesTable.(m){1}{:,:};

        frequencies = values(:,1)/1000;
        kz          = values(:,2);

        colourName = '';
        try
            colourName = colours{n+1,m};
        catch
        end
        
        if ~isempty(colourName)
            try
                colour = getColour(colourName);
            catch
                warning('Colour name %s not defined', colourName);
                colour = getColourByIndex(m);
            end   
        else
            colour = getColourByIndex(m);
        end
        
        plot(frequencies,kz, 'LineWidth', plotSettings.lines.thick, 'Color', colour);
        
    end
    hold off

    set(gca,'fontsize', plotSettings.labels.tickFontSize);
    
    xlim(fLim/1000)
    ylim(kLim)
    
    ylabel('$k_z~[rad/m]$',     'fontsize', plotSettings.labels.fontSize, 'interpreter', 'latex');
    xlabel('$Frequency~[kHz]$', 'fontsize', plotSettings.labels.fontSize, 'interpreter', 'latex');

    
    %legend
    if evanescent == true
        location = 'NorthEast';
    else
        location = 'NorthWest';
    end
    
    [legh,objh] = legend(lgd, 'Location',location,'FontSize', plotSettings.legend.fontSize);
    
    legh.BoxFace.ColorType  = plotSettings.legend.boxColourType;
    legh.BoxFace.ColorData  = plotSettings.legend.boxColourData;
    lineh                   = findobj(objh,'type','line');
    set(lineh,'LineWidth', plotSettings.legend.boxLineThickness);
     
    
    %figure
    set(modesFigure, 'Units', 'centimeters', 'Position', plotSettings.position, 'PaperUnits', 'centimeters', 'PaperSize', plotSettings.paperSize); 
    
 
    if overwrite
    
        if evanescent
            modeId =  [num2str(n), 'Xe',];
        else
            modeId =  [num2str(n), 'Xp'];
        end
            
        figureName = [ns.acMod.extractor.figureRootName, '_', modeId];
        saveas(modesFigure, [outputDirectory, '\', figureName, '.pdf'])

        close (modesFigure)
        
    end

end




% Save extracted data
%
% Inputs:  
%           modeRootDirectory: output folder
%           mModesTable: mode result table
%           n: mode order
%           evanescent: mode type (true/false)
%           overwrite: save and overwrite data (true/false)
%

function saveData(modeRootDirectory, nModesTable, n, evanescent, overwrite)

    if overwrite
        
        nameSpace      

        dataOutputFile      = [modeRootDirectory, '\', ns.acMod.extractor.resultsFileName, '.mat'];

        existingNmodeTable  = [];        
        if exist(dataOutputFile, 'file')
            existingNmodeTable = load(dataOutputFile);
            existingNmodeTable = existingNmodeTable.nModesTable;
        end

        if ~isempty(existingNmodeTable)

            %overwrite only the selected mode

            if evanescent
                modeTypeString = ns.acMod.modes.evanescentFolderName;
            else
                modeTypeString = ns.acMod.modes.propagativeFolderName;
            end

            if ~(sum(cellfun(@(x) strcmp(sprintf('n=%d',n), x), existingNmodeTable.Properties.VariableNames))>0)
            
                addCol             = nModesTable(:, strcmp(nModesTable.Properties.VariableNames, sprintf('n=%d',n)));
                existingNmodeTable = [existingNmodeTable, addCol];
                existingNmodeTable = existingNmodeTable(:,sort(existingNmodeTable.Properties.VariableNames));
                
            else
            
                existingNmodeTable.(sprintf('n=%d',n)).(modeTypeString) = nModesTable.(sprintf('n=%d',n)).(modeTypeString);

            end
            
            nModesTable = existingNmodeTable;
            
        end
    
        save(dataOutputFile, ns.acMod.extractor.outputDataFileName);
        
    end
        
end


        
    