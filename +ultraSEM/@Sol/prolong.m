function out = prolong(sol, n)
%PROLONG   Prolong an ULTRASEM.SOL.
%   PROLONG(SOL, N) returns an ULTRASEM.SOL representing the same function
%   as SOL, but using an N x N discretization on each element.

if ( nargin == 1 )
    out = sol;
    return
end

out = sol;
for k = 1:length(sol)
    [nx, ny] = size(sol.u{k});
    out.u{k} = zeros(n);
    ny = min(ny, n);
    nx = min(nx, n);
    out.u{k}(1:ny,1:nx) = sol.u{k}(1:ny,1:nx);
end

end
