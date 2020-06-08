function m = minEst(sol)
%MINEST   Estimate the minimum value of an ULTRASEM.SOL.
%
%   See also MAXEST, NORM.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

m = min(cellfun(@(u) min(u(:)), coeffs2vals(sol.coeffs)));

end
