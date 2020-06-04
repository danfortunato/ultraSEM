function m = minEst(sol)
%MINEST   Estimate the minimum value of an ULTRASEM.SOL.

m = min(cellfun(@(u) min(u(:)), coeffs2vals(sol.coeffs)));

end
