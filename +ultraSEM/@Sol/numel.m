function N = numel(sol)
%NUMEL   Number of degrees of freedom in an ULTRASEM.SOL.
%   N = NUMEL(S) returns the total number of degrees of freedom in the
%   ULTRASEM.SOL object S.

N = 0;
for k = 1:length(sol)
    N = N + numel(sol.coeffs{k});
end

end
