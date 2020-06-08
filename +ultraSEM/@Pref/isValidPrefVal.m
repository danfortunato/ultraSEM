function isValid = isValidPrefVal(prefName, prefVal)
%ISVALIDPREFVAL   Determine if a preference setting is valid.
%   ISVALID = ISVALIDPREFVAL(PREFNAME, PREFVAL) returns true if the value
%   PREFVALUE is a valid setting for the preference PREFNAME, and false
%   otherwise.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

stringPrefs = {'solver', 'discretization'};
functionPrefs = {'interfaceDegree'};
logicalPrefs = {'splitTriangles'};
integerPrefs = {'discSize'};

prefVals = ultraSEM.Pref.getValidPrefVals(prefName);

switch prefName
    case stringPrefs
        isValid = any(strcmp(prefVal, prefVals));
        if ( strcmp({prefName, prefVal}, {'discretization', 'values'}) )
            error('Discretization with values is not currently supported.');
        end
    case functionPrefs
        isValid = isa(prefVal, 'function_handle');
    case logicalPrefs
        isValid = islogical(prefVal) || any(prefVal == [0 1]);
    case integerPrefs
        isValid = isscalar(prefVal) && prefVal > 0 && prefVal == round(prefVal);
    otherwise
        error('ULTRASEM:PREF:isValidPrefVal:unknownPref', ...
                'Unknown preference name.');
end

end
