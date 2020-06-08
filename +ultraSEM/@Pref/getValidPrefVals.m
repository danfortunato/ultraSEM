function prefVals = getValidPrefVals(prefName)
%GETVALIDPREFVALS   Get valid preference values for a specific preference.
%   PREFVALS = ULTRASEM.PREF.GETVALIDPREFVALS(PREFNAME) returns a cell
%   array containing all known valid values of the preference PREFNAME.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

switch prefName
    case 'solver'
        prefVals = {'\', 'woodbury', 'LU'};
    case 'interfaceDegree'
        prefVals = {@max, @min, @mean};
    case 'splitTriangles'
        prefVals = {true, false};
    case 'discSize'
        preVals = {21};
    case 'discretization'
        prefVals = {'coeffs', 'values'};
    otherwise
        error('ULTRASEM:PREF:getValidPrefVals:unknownPref', ...
                'Unknown preference name.');
end

end
