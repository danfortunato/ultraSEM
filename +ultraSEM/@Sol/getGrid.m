function [x, y] = getGrid(sol, kk)
%GETGRID   Get the grid of an ULTRASEM.SOL.
%   [X, Y] = GETGRID(SOL) returns the tensor product Chebyshev grids (X, Y)
%   for the patches of SOL.
%
%   [X, Y] = GETGRID(SOL, KK) returns the tensor product Chebyshev grids
%   for the KK-th patches of SOL.

d = sol.domain;
u = sol.u;

x = cell(size(u));
y = cell(size(u));

if ( nargin == 1 )
    kk = 1:size(d, 1);
end

for k = kk
    nk = size(u{k}, 1);
    if ( isnumeric(d(k,:)) )
        [x{k,1}, y{k,1}] = util.chebpts2(nk, nk, d(k,:));
    else
        [xk, yk] = chebpts2(nk);
        [x{k,1}, y{k,1}] = transformGrid(d(k,:), xk, yk);
    end
end

if ( nargin > 1 && numel(kk) == 1 )
    x = x{1};
    y = y{1};
end

end
