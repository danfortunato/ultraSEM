function out = isInitialized(S)
%ISINITIALIZED   Check to see if an ULTRASEM has been initialized.

out = ~isempty(S.patches) && size(S.patches{1}.S, 2) > 0;

end
