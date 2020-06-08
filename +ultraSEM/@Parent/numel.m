function N = numel(P)
%NUMEL   Number of degrees of freedom in an ULTRASEM.PARENT.
%   N = NUMEL(P) returns the total number of degrees of freedom in the
%   ULTRASEM.PARENT object P.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

N = numel(P.child1) + numel(P.child2);

end
