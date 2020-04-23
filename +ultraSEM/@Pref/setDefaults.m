function setDefaults(varargin)
%SETDEFAULTS   Set default preferences.
%   ULTRASEM.PREF.SETDEFAULTS(PREF1, VAL1, PREF2, VAL2, ...) sets the
%   default values for the preferences whose names are stored in the
%   strings PREF1, PREF2, ..., etc. to VAL1, VAL2, ..., etc. All
%   subsequently constructed ULTRASEM.PREF objects will use these values
%   as the defaults.
%
%   ULTRASEM.PREF.SETDEFAULTS(PREF) sets the default values to the
%   preferences stored in the ULTRASEM.PREF object PREF. PREF can also be a
%   MATLAB structure, in which case it is converted to an ULTRASEM.PREF as
%   described in the documentation for the ULTRASEM.PREF constructor first.
%
%   ULTRASEM.PREF.SETDEFAULTS('factory') resets the default preferences to
%   their factory values.
%
% See also GETFACTORYDEFAULTS.

if ( nargin == 0 )
    ultraSEM.Pref.setDefaults('factory');
    return
end

if ( nargin == 1 )
    if ( isstruct(varargin{1}) )
        varargin{1} = ultraSEM.Pref(varargin{1});
    end

    if ( ischar(varargin{1}) && strcmp(varargin{1}, 'factory') )
        ultraSEM.Pref.manageDefaultPrefs('set-factory');
    elseif ( isa(varargin{1}, 'ultraSEM.Pref') )
        ultraSEM.Pref.manageDefaultPrefs('set', varargin{1}.prefList);
    else
        error('ULTRASEM:PREF:setDefaults:badArg', ...
            ['When calling ultraSEM.Pref.setDefaults() with just one ' ...
             'argument, that argument must be ''factory'' or a MATLAB ' ...
             'structure.']);
    end
elseif ( mod(nargin, 2) == 0 )
    ultraSEM.Pref.manageDefaultPrefs('set', varargin{:});
else
    error('ULTRASEM:PREF:setDefaults:unpairedArg', ...
        'Unpaired argument in name-value pair list.');
end

end
