function m = maxEst(sol)
%MAXEST   Estimate the maximum value of an ULTRASEM.SOL.

m = max(cellfun(@(u) max(u(:)), coeffs2vals(sol.coeffs)));

end
