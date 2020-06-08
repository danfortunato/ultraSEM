function pref = factoryDefaultPrefs()
%FACTORYDEFAULTPREFS   Get structure of factory default preferences.
%   S = ULTRASEM.PREF.FACTORYDEFAULTPREFS() returns a structure suitable
%   for storing in the prefList property of an ULTRASEM.PREF object that
%   contains all of the factory default values of the ULTRASEM preferences.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

pref.solver = 'woodbury';
pref.interfaceDegree = @max;
pref.splitTriangles = true;
pref.discSize = 21;
pref.discretization = 'coeffs';

end
