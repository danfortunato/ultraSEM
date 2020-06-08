function out = isInitialized(S)
%ISINITIALIZED   Check to see if an ULTRASEM has been initialized.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

out = ~isempty(S.patches) && size(S.patches{1}.S, 2) > 0;

end
