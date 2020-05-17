function out = isBuilt(S)
%ISBUILT   Check to see if an ULTRASEM has been built.

out = numel(S.patches) == 1;

end
