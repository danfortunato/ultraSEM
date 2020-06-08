function N = numel(L)
%NUMEL   Number of degrees of freedom in an ULTRASEM.LEAF.
%   N = NUMEL(L) returns the total number of degrees of freedom in the
%   ULTRASEM.LEAF object L.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

N = L.p^2;

end
