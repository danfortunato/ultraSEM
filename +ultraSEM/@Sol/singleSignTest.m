function [out, wzero] = singleSignTest(f)
%SINGLESIGNTEST   Heuristic check to see if F changes sign.
%   SINGLESIGNTEST(F) returns TRUE if the sampled values of the
%   ULTRASEM.SOL F have the same sign and FALSE otherwise.
%
%   [OUT, WZERO] = SINGLESIGNTEST(F) also returns a flag WZERO indicating
%   if any of the values are exactly zero.

% Evaluate on a grid:
X = cell2mat(coeffs2vals(f.coeffs));
X = X(:);

tol = eps;
vscale = norm(f, inf);
out = all( X >= -tol * vscale ) || ... % All values are nonnegative
      all( X <=  tol * vscale );       % All values are nonpositive

wzero = any(X == 0); % Any exact zeros on the grid?

end
