% SET VALUES IN IN THE TARGET STRUCT

% NOTES:
%
%   The fields of targetStruct are set with the values in settingsStruct.
%   If some field is missing, the values are retrieved from defaultStruct. 
%   Set ignoreMissing true to ingnore fields that are missing in both 
%   settingsStruct and defaultStruct.
%

% COPYRIGHT: 
%
%   MIT License
%   © 2023 Dario Chiantello <pipesoundlibrary@gmail.com>
%

function targetStruct = getSettings(targetStruct, settingsStruct, defaultStruct, ignoreMissing)

fieldNames      = fieldnames(targetStruct);
numberOfFields  = length(fieldNames);

for fieldCounter = 1:numberOfFields

    if isfield(settingsStruct, fieldNames{fieldCounter})
        
        targetStruct.(fieldNames{fieldCounter}) = settingsStruct.(fieldNames{fieldCounter});
               
    elseif isfield(defaultStruct, fieldNames{fieldCounter})
        
        %warning('Default value used for field %s', fieldNames{fieldCounter})
        targetStruct.(fieldNames{fieldCounter}) = defaultStruct.(fieldNames{fieldCounter});
        
    else
        if ~ignoreMissing
        
            error('Setting field %s not found', fieldNames{fieldCounter});
            
        end
    end
    
end
    
   