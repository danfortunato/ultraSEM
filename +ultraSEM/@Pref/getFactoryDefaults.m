function pref = getFactoryDefaults()
%GETFACTORYDEFAULTS   Get factory default preferences.
%   PREF = ULTRASEM.PREF.GETFACTORYDEFAULTS() returns an ULTRASEM.PREF
%   object with the preferences set to their factory defaults, irrespective
%   of the currently defined values of the default preferences. This
%   function is useful if the user wishes to solve PDEs with ULTRASEM using
%   the factory defaults when other user-set defaults are currently in
%   force.
%
% See also SETDEFAULTS.

fd = ultraSEM.Pref.factoryDefaultPrefs();
pref = ultraSEM.Pref(fd);

end
