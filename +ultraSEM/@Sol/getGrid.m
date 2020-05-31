function [x, y] = getGrid(sol, kk)
%GETGRID   Get the grid of an ULTRASEM.SOL.
%   [X, Y] = GETGRID(SOL) returns the tensor product Chebyshev grids (X, Y)
%   for the patches of SOL.
%
%   [X, Y] = GETGRID(SOL, KK) returns the tensor product Chebyshev grids
%   for the KK-th patches of SOL.

d = sol.domain;
u = sol.u;

if ( nargin == 1 )
    kk = 1:size(d, 1);
end

x = cell(numel(kk),1);
y = cell(numel(kk),1);

i = 1;
for k = kk
    nk = size(u{k}, 1);
    if ( isnumeric(d(k,:)) )
        [x{i,1}, y{i,1}] = util.chebpts2(nk, nk, d(k,:));
    else
        [xk, yk] = chebpts2(nk);
        [x{i,1}, y{i,1}] = transformGrid(d(k,:), xk, yk);
    end
    i = i+1;
end

if ( nargin > 1 && numel(kk) == 1 )
    x = x{1};
    y = y{1};
end

end
