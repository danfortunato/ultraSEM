function display(pref) %#ok<DISPLAY>
%DISPLAY   Display an ULTRASEM.PREF object.
%   DISPLAY(PREF) prints out a list of the preferences stored in the
%   ULTRASEM.PREF object PREF.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

% Compute the screen column in which pref values start.
col = 20;
ind = '    ';

% Print values of "known" preferences.
prefList = pref.prefList;

fprintf('ultraSEM.Pref object with the following preferences:\n');
fprintf([ind pad('discretization:', col) '%s\n'], ...
    prefList.discretization);
fprintf([ind pad('discSize:', col) '%d\n'], ...
    prefList.discSize);
fprintf([ind pad('interfaceDegree:', col) '%s\n'], ...
    func2str(prefList.interfaceDegree));
fprintf([ind pad('solver:', col) '%s\n'], ...
    prefList.solver);
fprintf([ind pad('splitTriangles:', col) '%s\n'], ...
    logical2str(prefList.splitTriangles));

end

function s = logical2str(x)
%LOGICAL2STR   Convert a logical to 'true' or 'false'.
    str = {'false', 'true'};
    s = str{x+1};
end
