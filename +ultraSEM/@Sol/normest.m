function val = normest(sol)
%NORMEST   Estimate maximum absolute value of SOL over all patches.

val = max(cellfun(@(u) max(abs(u(:))), coeffs2vals(sol.u)));

end
