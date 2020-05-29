function [x, y] = plotpts(sol, kk)
%PLOTPTS   Get plot points.
%   [X, Y] = PLOTPTS(SOL) returns points (X, Y) which are uniformly spaced
%   over each patch in SOL, to be used for plotting.
%
%   [X, Y] = PLOTPTS(SOL, KK) returns plot points for the KK-th patches.

d = sol.domain;
coeffs = sol.coeffs;

x = cell(size(coeffs));
y = cell(size(coeffs));

if ( nargin < 2 )
    kk = 1:size(d, 1);
end

kk = kk(:).';

for k = kk
    if ( isnumeric(d(k,:)) )
        [x{k,1}, y{k,1}] = meshgrid(linspace(d(k,1), d(k,2), sol.nplotpts), ...
                                    linspace(d(k,3), d(k,4), sol.nplotpts));
    else
        [xk, yk] = meshgrid(linspace(-1, 1, sol.nplotpts));
        [x{k,1}, y{k,1}] = transformGrid(d(k,:), xk, yk);
    end
end

if ( nargin > 1 && numel(kk) == 1 )
    x = x{1};
    y = y{1};
end

end
