function out = isempty(sol)
%ISEMPTY   Test for empty ULTRASEM.SOL.
%   ISEMPTY(SOL) returns 1 if SOL is an empty ULTRASEM.SOL and 0 otherwise.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

out = isempty(sol.coeffs);

end
