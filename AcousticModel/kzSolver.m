% NUMERICAL SOLVER FOR THE EXTRACTION OF THE AXIAL WAVENUMBERS kz

% COPYRIGHT: 
%
%   MIT License
%   © 2023 Dario Chiantello <pipesoundlibrary@gmail.com>
%

% DESCRIPTION:
% 
%   Numerical solver for the calculation of the axial wavenumber kz for an
%   acoustic waveguide composed of a cylindrical elastic shield filled with
%   inviscid liquid. Modes can be either propagative or evanescent, but no 
%   attenuation is considered. The solver provides the solution within a
%   predefined frequency range, for a given kz range, and for a given
%   angular mode order. The solver generates an intermediate results file 
%   and a graphic representation of the solution.
%

% INPUTS: 
%   
%   Input settings can be provided through external setting files or
%   obtained from spaceSettings.m as default values. The solver also 
%   requires the AcusticModel/characteristicMatrix.m, which defines the
%   characteristic matrix of the waveguide. The customised settings 
%   definition requires two files stored under 
%   [simulation folder]\SimulationSetup:
% 
%   - simulationSetup.m defines the geometry and the materials of
%     the simulation,
%   - solverSetup.m defines the solver settings. 
% 
%   Examples of setting files can be found in the examples for the acoustic
%   model. The setting structs are reported below, along with the related 
%   description. Filed not specified are obtained from spaceSettings.m as 
%   default values. 
%
%
%   simulationSetup.m
%
%   -sis.
%      R1_m                             Internal radius of the pipe in m                 
%      R2_m                         	External radius of the pipe in m           
%      c_phi_m_s                        Sound velocity of the inner liquid in m/s                
%      ro_l_kg_m3                       Density of the inner liquid in kg/m3              
%      c_gamma_m_s                      Longitudinal velocity of the solid shield in m/s                   
%      c_psi_m_s                        Shear velocity of the solid shield in m/s                  
%      ro_s_kg_m3                       Density of the solid shield in kg/m3                
%
%
%   solverSetup.m
%
%   -sos.settings
%      nMin                             Min angular mode order n
%      nMax                             Max angular mode order n (N)
%      evanProp                         0-only propagative, 1-only evanescent, 2-both
%      kMin_rad_m                       kz min range rad/m               
%      kMax_rad_m                       kz max range rad/m
%      kResolution_rad_m                kz resolution rad/m   
%      kSearchRadius_rad_m              kz serach radius rad/m
%      kDiscontinuityRadius_rad_m       kz radius rad/m for the discontinuity check
%      cutFrequencySearchScale          Finer scale of frequency resolution for the determination of the kz=0 cut frequencies
%      fDiscontinuityRadius_Hz          f discontinuity radius in Hz for the discontinuity check of the cut-off
%      fMin_Hz                          Frequency min range in Hz             
%      fMax_Hz                          Frequency max range Hz             
%      fResolution_Hz                   Frequency resolution Hz   
%      fInitSteps                       Number of initial frequency steps
%      fZeroNeighboursSteps             Number of neighbouring points to check in order to determine the kz search intervals
%      zeroCutExtraInitFreqLeap_Hz      Frequency distance in Hz between extra init points around the kz=0 cut frequencies
%      zeroCutExtraInitPointsRadius     Number of extra init points (per side) around kz=0 cut frequencies
%      fExtraInitSteps_Hz               Array of extra initialisation frequencies in Hz (e.g. [1035, 1498])
% 
%   -sos.pictures.
%      format                           Picture format, see Utilities/getPlotSettings 
%      size                             Picture size, see Utilities/getPlotSettings 
%

% OUTPUTS: 
% 
%   The outputs of the solver are the kz values for the waveguide defined. 
%   The following files are generated:
%   - solverResult_XXy.mat: Intermediate processing results generated by 
%     the solver. Placed in [simulation folder]\Modes\Mode_XX\modeType.
%   - solverResult_XXy.pdf: Graphical representation of the raw solution.
%     Placed in [simulation folder]\Modes\Mode_XX\modeType
% 

% PROCESSING STEPS: 
% 
%   1- The first step yields the zero-cut frequencies obtained for kz=0.
% 
%   2- The second step initialises the search by scanning along the defined 
%      kz spectrum. This scan is executed for a limited number of 
%      frequencies at regular intervals and around the zero-cut frequencies 
%      found in the previous step. 
% 
%   3- The third step scans along the whole frequency spectrum, focusing 
%      only in the neighbourhood of the kz points already found. 
% 

% NOTES: 
% 
%   No settings to edit in this file. Change the name of the simulation 
%   folder if required. The root repository is defined by the field
%   ns.repository.path in nameSpace.m. To change the default settings, add 
%   the customised setting files under [simulationFolder]\SimulationSetup
%


clear
clc


%DEFAULT SIMULATION NAME
defaultSimulationName = 'ALU R1 0984 R2 1000'; %test%;


%LOAD LIBRARY
loadLibray();
nameSpace
spaceSettings


%NAME OF THE SIMULATION FOLDER
prompt          = sprintf("Default simulation name: ""%s"". Press enter to confirm or type a new name.\n", defaultSimulationName);
simulationName  = input(prompt, "s");

if isempty(simulationName)
    
    simulationFolderName = defaultSimulationName;
    
else
    
    simulationFolderName = simulationName;
    
end


%SAVE/OVERWRITE DATA
overwrite      = false;
prompt         = sprintf("Do you want to save and overwrite solver result data? [Y/N] >\n");
overwriteAns   = input(prompt, "s");
if strcmp(overwriteAns, 'Y') ||  strcmp(overwriteAns, 'y')
    overwrite = true;
else
    warning("New data won't not be saved");
end



%SIMULATION FOLDER  
simulationFolderPath                    = [ns.repository.path, '\', ns.acMod.folderName,'\',simulationFolderName];
if ~exist(simulationFolderPath, 'dir')
   mkdir(simulationFolderPath)
end


%MODE ROOT FOLDER
modeRootDirectory                       = [simulationFolderPath,'\', ns.acMod.modes.folderName];
if ~exist(modeRootDirectory, 'dir')
   mkdir(modeRootDirectory)
end


%SIMULATION AND CALCULATION SETUP           
simulationSetupFolderPath               = [simulationFolderPath,'\',ns.acMod.setup.folderName];
if ~exist(simulationSetupFolderPath, 'dir')
   mkdir(simulationSetupFolderPath)
end

addpath(simulationSetupFolderPath);

if exist([simulationSetupFolderPath, '\', ns.acMod.setup.simulation.fileName, '.m'], 'file')
    
    sisIn = simulationSetup();
    
else
    
    wString = sprintf('Simulation setup file missing in: %s\n', simulationSetupFolderPath);
    wString = [wString, 'The solver will use the default values defined in spaceSettins.m'];
    warning(wString)    
    sisIn = [];
    
end

if exist([simulationSetupFolderPath, '\', ns.acMod.setup.kzSolver.fileName, '.m'], 'file')
    
    sosIn = solverSetup(0,0);
   
else
    
    wString = sprintf('Solver setup file missing in: %s\n', simulationSetupFolderPath);  
    wString = [wString, 'The solver will use the default values defined in spaceSettins.m'];
    warning(wString) 
    sosIn = [];
    
end

%simulation settings
sis.R1_m                                        = [];
sis.R2_m                                        = [];
sis.c_phi_m_s                                   = [];
sis.ro_l_kg_m3                                  = [];
sis.c_gamma_m_s                                 = [];   
sis.c_psi_m_s                                   = [];
sis.ro_s_kg_m3                                  = [];

if ~isfield(ss.acMod, 'simSettings')  
    ss.acMod.simSettings = [];
end
sis = getSettings(sis, sisIn, ss.acMod.simSettings, false);


%solver setup
sos.settings.nMin                                = [];                                         
sos.settings.nMax                                = [];
sos.settings.evanProp                            = [];
sos.settings.kMin_rad_m                          = [];
sos.settings.kMax_rad_m                          = [];
sos.settings.kResolution_rad_m                   = [];
sos.settings.kSearchRadius_rad_m                 = [];
sos.settings.kDiscontinuityRadius_rad_m          = [];
sos.settings.cutFrequencySearchScale             = [];
sos.settings.fDiscontinuityRadius_Hz             = [];
sos.settings.fMin_Hz                             = [];
sos.settings.fMax_Hz                             = [];
sos.settings.fResolution_Hz                      = [];
sos.settings.fInitSteps                          = [];
sos.settings.fZeroNeighboursSteps                = [];
sos.settings.zeroCutExtraInitFreqLeap_Hz         = [];
sos.settings.zeroCutExtraInitPointsRadius        = [];
sos.settings.fExtraInitSteps_Hz                  = [];

if ~isfield(ss.acMod, 'solver')
    ss.acMod.solver = [];
end
if ~isfield(sosIn, 'settings')
    sosIn.settings = [];
end
sos.settings = getSettings(sos.settings, sosIn.settings, ss.acMod.solver, false);
    

%figure setup
sos.pictures.format = '';
sos.pictures.size   = '';
if isfield(sosIn, 'plot')

    if isfield(sosIn.pictures, 'format')
       sos.pictures.format = sosIn.pictures.format;
    end
        
    if isfield(sosIn.pictures, 'size')
       sos.pictures.size = sosIn.pictures.size;
    end
end
    



%SOLVER ROUTINE
if sos.settings.evanProp <= 1
    evanProp = sos.settings.evanProp;
else
    evanProp = 0:1;
end

%LOOP OVER THE MODES
for n = sos.settings.nMin:sos.settings.nMax
        
    for evanescent = evanProp
        
        sosIn        = solverSetup(n, evanescent);
        if ~isfield(sosIn, 'settings')
            sosIn.settings = [];
        end
        sos.settings = getSettings(sos.settings, sosIn.settings, ss.acMod.solver, true);
        
        modeDirectoryName = [ns.acMod.modes.modeFolderRootName, '_', num2str(n), 'X'];

        if evanescent == true   
            modeType = ns.acMod.modes.evanescentFolderName; 
        else
            modeType = ns.acMod.modes.propagativeFolderName;  
        end

        outputDirectory = [modeRootDirectory, '\', modeDirectoryName, '\', modeType];
        if ~exist(outputDirectory, 'dir')
            mkdir(outputDirectory);
        end

        
        tic()
        fprintf('\n**********\n')        
        fprintf("Calculating n=%d %s modes...\n", n, modeType)

        
        %ZERO CUT FREQUENCIES
        fprintf("Calculating zero cut frequencies...\n")
        
        kz                           = 0;
        fCutFrequenciesVector        = sos.settings.fMin_Hz : round(sos.settings.fResolution_Hz/sos.settings.cutFrequencySearchScale) : sos.settings.fMax_Hz;
        
        zeroKcutFrequencyComplex     = zeros(length(fCutFrequenciesVector),1); 
        omega                        = 2*pi* fCutFrequenciesVector;
        
        for fCutFrequencyCounter = 1:length(fCutFrequenciesVector)     

            zeroKcutFrequencyComplex(fCutFrequencyCounter) = det(characteristicMatrix(n, omega(fCutFrequencyCounter), kz, sis)); 
            
        end
        
        zeroKcutFrequencyReal        = real(zeroKcutFrequencyComplex);
        zeroKcutFrequencyRealCheck   = zeroKcutFrequencyReal > 0;
        zeroKcutFrequencyRealCheck   = xor(zeroKcutFrequencyRealCheck, circshift(zeroKcutFrequencyRealCheck,1));

        cutFreqDiscontinuityRadius   = round(sos.settings.fDiscontinuityRadius_Hz/(sos.settings.fResolution_Hz/sos.settings.cutFrequencySearchScale));
        zeroCheckCutFrequencyVector  = discontinuityCheck(zeroKcutFrequencyRealCheck, abs(zeroKcutFrequencyReal), cutFreqDiscontinuityRadius);

        %save cut frequency table
        zeroCutFrequencies           = fCutFrequenciesVector(zeroCheckCutFrequencyVector == 1);        
        
        fprintf("Done!\n")
        
        
        
        %Kz
        fprintf("Calculating kz...\n")
        
        
        %init vars        
        kVector                     = sos.settings.kMin_rad_m : sos.settings.kResolution_rad_m : sos.settings.kMax_rad_m;
        fVector                     = sos.settings.fMin_Hz : sos.settings.fResolution_Hz : sos.settings.fMax_Hz;
        fInitLeap                   = round(length(fVector)/sos.settings.fInitSteps);
        fInitVectorIds              = 1 : fInitLeap : length(fVector);
        discontinuityCheckPoints    = round(sos.settings.kDiscontinuityRadius_rad_m/sos.settings.kResolution_rad_m);
        zeroCheckMatrix             = false(length(kVector), length(fVector), 2);
        omega                       = 2*pi* fVector;
        
        if evanescent == true
            ikVector = i*kVector;
        else
            ikVector = kVector;   
        end 
        
        
        %add extra init steps (user defined and around the kz=0 cut frequencies )
        [~,extraIds]    = arrayfun(@(x) min(abs(fVector - x)), sos.settings.fExtraInitSteps_Hz, 'UniformOutput', false);
        extraIds        = cell2mat(extraIds);
        
        zeroCutExtraIds = [];
        for zeroCutCounter = 1:length(zeroCutFrequencies)
            
            zeroCutExtraInitFreq = zeroCutFrequencies(zeroCutCounter) - sos.settings.zeroCutExtraInitPointsRadius*sos.settings.zeroCutExtraInitFreqLeap_Hz :...
                                   sos.settings.zeroCutExtraInitFreqLeap_Hz:...
                                   zeroCutFrequencies(zeroCutCounter) + sos.settings.zeroCutExtraInitPointsRadius*sos.settings.zeroCutExtraInitFreqLeap_Hz;
            
           [~,extraPoints]  = arrayfun(@(x) min(abs(fVector - x)), zeroCutExtraInitFreq, 'UniformOutput', false);
            
            zeroCutExtraIds = [zeroCutExtraIds, extraPoints];
           
        end
        zeroCutExtraIds = cell2mat(zeroCutExtraIds);
        

        fInitVectorIds  = [fInitVectorIds, extraIds, zeroCutExtraIds];
        fInitVectorIds  = unique(fInitVectorIds);
        fInitVectorIds  = sort([fInitVectorIds, extraIds]);

        
        fInitVector     = fVector(fInitVectorIds);
        
        leftSteps       = fInitVectorIds - circshift(fInitVectorIds, 1) - 1;
        leftSteps(1)    = fInitVectorIds(1) - 1;
        rightSteps      = circshift(fInitVectorIds, -1) - fInitVectorIds - 1;
        rightSteps(end) = length(fVector) - fInitVectorIds(end);
        scanSteps       = zeros(2, length(fInitVector));    
        scanSteps(1,:)  = leftSteps;
        scanSteps(2,:)  = rightSteps;
        

        completionCounter   = 0;
        percentageStringLen = fprintf('Percentage completed: %d%%\n', round(completionCounter));
        
        %scan points starting from the initialisation points
        for fInitCounter = 1:length(fInitVectorIds)     
        
            %update completion indicator
            completionCounter = fInitCounter/length(fInitVector) * 100;
            fprintf(repmat('\b',1,percentageStringLen));             
            percentageStringLen = fprintf('Percentage completed: %d%%\n', round(completionCounter));
            
            
            %init at initialisation point
            fVectorPointer  = fInitVectorIds(fInitCounter);
                    
            zeroComplex     = zeros(length(kVector), 1);
            
            for kzInitCounter = 1:length(kVector)
            
                zeroComplex(kzInitCounter) = det(characteristicMatrix(n, omega(fVectorPointer), ikVector(kzInitCounter), sis));
                
            end
            
            zeroReal                                = real(zeroComplex);
            zeroRealCheck                           = zeroReal > 0;     
            zeroRealCheck                           = xor(zeroRealCheck, circshift(zeroRealCheck,1));
            zeroRealCheck                           = discontinuityCheck(zeroRealCheck, abs(zeroReal), discontinuityCheckPoints);
            zeroRealCheck(1)                        = 0;
            zeroRealCheck(end)                      = 0;   
            zeroCheckMatrix(:, fVectorPointer,1)    = zeroRealCheck;
            zeroCheckMatrix(:, fVectorPointer,2)    = zeroRealCheck;
            
            
            %scan at both sides of the initialisation point (1-left 2-right)
            for scanDirection = 1:2 
                                
                fVectorInitPointer = fInitVectorIds(fInitCounter);
                
                for scanCounter = 1: scanSteps(scanDirection, fInitCounter)

                    if scanDirection == 1

                        fVectorPointer      = fVectorInitPointer - scanCounter;
                        
                        maxNeighbour        = min(fVectorPointer+sos.settings.fZeroNeighboursSteps, length(fVector));
                        neighbourhood       = fVectorPointer+1 : maxNeighbour;
                        
                    else

                        fVectorPointer      = fVectorInitPointer + scanCounter;
                        
                        minNeighbour        = max(fVectorPointer-sos.settings.fZeroNeighboursSteps, 1);
                        neighbourhood       = minNeighbour : fVectorPointer-1;
                        
                    end

                    %get the neighbouring points where to search
                    zeroPointsIds = [];
                    for neighbour = neighbourhood    

                        zeroPointsIds = [zeroPointsIds, find(zeroCheckMatrix(:, neighbour, scanDirection))'];

                    end
                    zeroPointsIds = unique(zeroPointsIds)';
                    
                    
                    topKsearch          = 0;
                    checkedValuesPos    = zeros(length(kVector),1);

                    zeroComplex = zeros(length(kVector), 1);
                    
                    for zeroPointCounter = 1:length(zeroPointsIds)

                        searchRadius    = round(sos.settings.kSearchRadius_rad_m/sos.settings.kResolution_rad_m);
                        botKsearch      = max(1, zeroPointsIds(zeroPointCounter) - searchRadius);
                        botKsearch      = max(botKsearch, topKsearch+1);
                        topKsearch      = min(length(kVector), zeroPointsIds(zeroPointCounter) + searchRadius);
                        
                        checkedValuesPos(botKsearch:topKsearch) = 1;
                                                
                        for kSearchId = botKsearch:topKsearch

                            zeroComplex(kSearchId) = det(characteristicMatrix(n, omega(fVectorPointer), ikVector(kSearchId), sis));

                        end

                    end

                    %remove the boundaries of the search region
                    botBoundaries = find ( (checkedValuesPos - circshift(checkedValuesPos, 1)) == 1);
                    topBoundaries = find ( (checkedValuesPos - circshift(checkedValuesPos,-1)) == 1);
                    checkedValuesPos(botBoundaries) = 0;
                    checkedValuesPos(topBoundaries) = 0;
                    
                    
                    zeroReal     	= real(zeroComplex);
                    zeroRealCheck   = zeroReal > 0;
                    zeroRealCheck   = xor(zeroRealCheck, circshift(zeroRealCheck,1));
                    zeroRealCheck(~checkedValuesPos) = 0;
                    zeroRealCheck   = discontinuityCheck(zeroRealCheck, abs(zeroReal), discontinuityCheckPoints);
                    zeroCheckMatrix(:, fVectorPointer, scanDirection) = zeroRealCheck;

                end
            end
        end
            
        %merge both scan directions
        zeroCheckMatrix = or(zeroCheckMatrix(:,:,1), zeroCheckMatrix(:,:,2));
        
        
        fprintf("Done!\n")
        toc();
        
        
        
        %OUTPUTS
        
        %save data        
        if evanescent
            modeId =  [num2str(n), 'Xe',];
        else
            modeId =  [num2str(n), 'Xp'];
        end
        
        outputFileName  = [ns.acMod.kzSolver.resultRootName, '_', modeId];
        outputFile      = [outputDirectory, '\', outputFileName];
        
        compactModeMatrix = compressMatrix(zeroCheckMatrix, fVector, kVector);
        
        if overwrite
            
            save(   outputFile,...
                    'compactModeMatrix',...
                    'zeroCutFrequencies',...
                    'n',...
                    'evanescent',...
                    'sis',... 
                    'sos',...
                    'fVector',...
                    'kVector')
        end
                
        
        %print and save figure
        printFigure(outputDirectory, zeroCheckMatrix, compactModeMatrix, n, evanescent, fVector, kVector, fInitVector, sos.pictures, overwrite)
        
   end           
end   




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



% Check for change of sign given by discontinuities and remove false zeros in zeroCheckMatrix
%
% Inputs: 
%         zeroVector: bool col vector indication the position of the zeros;
%         modulusVector: col vector of abs(det(A));
%         discontinuityRadius: radius for the local maximum check
%
% Output: 
%         zeroVector: same as input but cleaned from false zeros

function zeroVector = discontinuityCheck(zeroVector, modulusVector, discontinuityRadius)

    foundZerosIds = find(zeroVector);

    for foundZerosCounter = 1:length(foundZerosIds)

        botttomRange                    = max(1,                    foundZerosIds(foundZerosCounter) - discontinuityRadius);
        topRange                        = min(length(zeroVector),   foundZerosIds(foundZerosCounter) + discontinuityRadius);
 
        localAbs                        = modulusVector;
        localAbs(1:botttomRange)        = 0;
        localAbs(topRange:end)          = 0;

        [~,localMaxId]                  = max(localAbs);

        %discontinuity found: in case of discontinuity the max is where the sign changes
        if( abs(localMaxId - foundZerosIds(foundZerosCounter)) <=2 )

            zeroVector(foundZerosIds(foundZerosCounter)) = 0;

        end
    end
end
    



%Compress the matrix before saving
%
% Input:  
%         fullMatrix: mode matrix to compress
%         fVector: frequency step vector
%         kVector: kz step vector
%
% Output: 
%         compactMatrix: compact mode matrix

function compactMatrix = compressMatrix(fullMatrix, fVector, kVector)

    compactMatrix           = zeros(sum(sum(fullMatrix)),2);
    compactMatrixCounter    = 0;
        
    for fCounter = 1 : length(fVector)
                    
        for kCounter = 1:length(kVector)

            if fullMatrix(kCounter, fCounter) == 1

                compactMatrixCounter = compactMatrixCounter+1;
                compactMatrix(compactMatrixCounter,:) = [fVector(fCounter),kVector(kCounter)];

            end
        end
    end
    
    kRes = 2*(kVector(2)-kVector(1));
    
    for fCounter = 1 : length(fVector)
    
        checkValuesIds = find(compactMatrix(:,1) == fVector(fCounter));
        checkValues    = compactMatrix(checkValuesIds,:);    
        
        for j = 1:length(checkValuesIds)
        
            deleteIds      = abs(checkValues(j,2) - checkValues(:,2)) < kRes;
            deleteIds(1:j) = 0;
            
            compactMatrix(checkValuesIds(deleteIds), :) = 0;
        end
    end
        
    compactMatrix = compactMatrix(compactMatrix(:,1)>0, :);
    
end




%Print and save the solver image
%
% Inputs: 
%           resultFolder: output result folder
%           resultMatrix: binary result matrix
%           compactModeMatrix: compact result matrix
%           n: mode order
%           evanescent: mode type (true/false)
%           fVector: frequency vector
%           kVector: kz test points vector
%           fInitVector: initialisation frequencies
%           pictureSettings: figure settings
%           overwrite: save and overwrite data (true/false)

function printFigure(resultFolder, resultMatrix, compactModeMatrix, n, evanescent, fVector, kVector, fInitVector, pictureSettings, overwrite)
     
    nameSpace

    plotSettings    = getPlotSettings(pictureSettings);
        
    thicknessX      = round(size(resultMatrix,2)/1000);
    thicknessY      = round(size(resultMatrix,1)/1000);
    imageBold       = createBoldImage(resultMatrix, thicknessX, thicknessY);
    
    paperSize       = plotSettings.paperSize;
    [h w]           = size(imageBold);
    w               = h * paperSize(1) / paperSize(2);

    maxImageLen     = 5000;

    if w > maxImageLen || h > maxImageLen

        scaleFactor = max(h/maxImageLen, w/maxImageLen);
        w = w/scaleFactor;
        h = h/scaleFactor;

    end

    resizedImage = imresize(imageBold, [h w], 'nearest');

    if evanescent
        figureNameString = sprintf('Evenescent modes n=%d\n', n);
    else
        figureNameString = sprintf('Propagative modes n=%d\n', n);
    end
    
    modesFigure = figure('Name', figureNameString);

    image(255*(~resizedImage))
    colormap(gray)
    set(gca,'YDir','normal');

    hold on
    axis on;
    box on

    xticks( 0:w/10:w )
    fStartKhz = fVector(1)/1000;
    fEndKhz   = fVector(end)/1000;
    xticklabels( arrayfun(@(x) sprintf('%.1f',x), fStartKhz:(fEndKhz-fStartKhz)/10:fEndKhz, 'UniformOutput', false ) );

    yticks( 0:h/8:h )
    yticklabels( arrayfun(@(x) sprintf('%.1f',x), kVector(1):(kVector(end)-kVector(1))/8:kVector(end), 'UniformOutput', false ) );

    set(gca,'fontsize', plotSettings.labels.tickFontSize);
    
    ylabel('$k_z~[rad/m]$',     'fontsize', plotSettings.labels.fontSize, 'interpreter','latex');
    xlabel('$Frequency~[kHz]$', 'fontsize', plotSettings.labels.fontSize, 'interpreter','latex');

    xlim([0 w])
    ylim([0 h])

    %red circles
    for pointCounter = 1:size(compactModeMatrix, 1)

        if sum(compactModeMatrix(pointCounter, 1) == fInitVector) > 0

            xPos = (compactModeMatrix(pointCounter,1) - fVector(1))/(fVector(end) - fVector(1))*w;
            yPos = (compactModeMatrix(pointCounter,2) - kVector(1))/(kVector(end) - kVector(1))*h;

            plot(xPos, yPos, 'o', 'LineWidth', getPlotLineThickness('thickLine'), 'Color', getColour('red'))

        end

    end
    
    %green lines
    for fLineCounter = 1:length(fInitVector)

        xline((fInitVector(fLineCounter) - fVector(1)) /(fVector(end) - fVector(1))*w, '--', 'LineWidth', getPlotLineThickness('thinLine'), 'Color', getColour('green') );

    end
    

    hold off

    set(modesFigure, 'Units', 'centimeters', 'Position', plotSettings.position, 'PaperUnits', 'centimeters', 'PaperSize', paperSize); 
 
    if overwrite
    
        if evanescent
            modeId =  [num2str(n), 'Xe',];
        else
            modeId =  [num2str(n), 'Xp'];
        end
        
        outputFileName  = [ns.acMod.kzSolver.resultRootName, '_', modeId];
        outputFile      = [resultFolder, '\', outputFileName, '.pdf'];
        
        saveas(modesFigure, outputFile)
        close(modesFigure)

    end
    
end

    


% Highlight the lines of the plot obtained from the binary matrix
%
% Inputs: 
%         imageMatrix: input image matrix;
%         thicknessX: thickness X;
%         thicknessY: thickness Y
%
% Output: 
%         imageBoldMatrix: output image matrix

function imageBoldMatrix = createBoldImage(imageMatrix, thicknessX, thicknessY)

    imageBoldMatrix = imageMatrix;
    temp            = imageMatrix;

    for i = -thicknessX:thicknessX

        imageBoldMatrix = or(imageBoldMatrix, circshift(temp,i,2));

    end

    temp = imageBoldMatrix;

    for i = -thicknessY:thicknessY

        imageBoldMatrix = or(imageBoldMatrix, circshift(temp,i,1));    

    end  
    
end
    
    
    
    
