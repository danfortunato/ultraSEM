function N = numel(L)
%NUMEL   Number of degrees of freedom in an ULTRASEM.LEAF.
%   N = NUMEL(L) returns the total number of degrees of freedom in the
%   ULTRASEM.LEAF object L.

N = L.p^2;

end
