function out = isBuilt(S)
%ISBUILT   Check to see if an ULTRASEM has been built.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

out = numel(S.patches) == 1;

end
