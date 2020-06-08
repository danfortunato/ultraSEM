function out = length(T)
%LENGTH   The length of an ULTRASEM.DOMAIN is its number of patches.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

out = size(T.domain, 1);

end
