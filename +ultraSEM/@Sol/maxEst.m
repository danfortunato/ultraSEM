function m = maxEst(sol)
%MAXEST   Estimate the maximum value of an ULTRASEM.SOL.
%
%   See also MINEST, NORM.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

m = max(cellfun(@(u) max(u(:)), coeffs2vals(sol.coeffs)));

end
