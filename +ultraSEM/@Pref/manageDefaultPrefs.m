function varargout = manageDefaultPrefs(varargin)
%MANAGEDEFAULTPREFS   Private method for handling default preferences.
%   ULTRASEM.PREF.MANAGEDEFAULTPREFS('get') returns a structure suitable
%   for storing in the prefList property of an ULTRASEM.PREF with all of
%   the currently stored default preferences suitable for initializing an
%   ULTRASEM.PREF object.
%
%   ULTRASEM.PREF.MANAGEDEFAULTPREFS('set-factory') restores the default
%   preferences to their factory values.
%
%   ULTRASEM.PREF.MANAGEDEFAULTPREFS('set', PREFLIST) sets the default
%   values to those stored in the structure PREFLIST. PREFLIST should be a
%   structure suitable for use as an ULTRASEM.PREF prefList.
%
%   ULTRASEM.PREF.MANAGEDEFAULTPREFS('set', PREF1, VAL1, PREF2, VAL2, ...)
%   sets the default values for PREF1, PREF2, ..., etc. to VAL1, VAL2, ...,
%   etc.

persistent defaultPrefs;

if ( isempty(defaultPrefs) )
    defaultPrefs = ultraSEM.Pref.factoryDefaultPrefs();
end

if ( strcmp(varargin{1}, 'get') )
    varargout{1} = defaultPrefs;
elseif ( strcmp(varargin{1}, 'set-factory') )
    defaultPrefs = ultraSEM.Pref.factoryDefaultPrefs();
elseif ( strcmp(varargin{1}, 'set') )
        varargin(1) = [];
    if ( isstruct(varargin{1}) )
        defaultPrefs = varargin{1};
    else
        while ( ~isempty(varargin) )
            prefName  = varargin{1};
            prefValue = varargin{2};
            if ( isfield(defaultPrefs, prefName) )
                if ( ultraSEM.Pref.isValidPrefVal(prefName, prefValue) )
                    defaultPrefs.(prefName) = prefValue;
                else
                    error('ULTRASEM:PREF:manageDefaultPrefs:invalidPrefVal', ...
                        'Invalid preference value.');
                end
            else
                error('ULTRASEM:PREF:manageDefaultPrefs:badPref', ...
                    'Unrecognized preference name.');
            end
            varargin(1:2) = [];
        end
    end
end

end
