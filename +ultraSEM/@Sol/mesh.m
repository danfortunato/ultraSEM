function varargout = mesh(sol, varargin)
%MESH   3-D mesh surface of an ULTRASEM.SOL.
%   MESH(SOL) plots the tensor product Chebyshev grids for the patches of
%   SOL, colored according to the values of SOL.
%
% See also PLOT, SURF, CONTOUR, GETGRID.

holdState = ishold();

vals = coeffs2vals(sol.u);
[x,y] = getGrid(sol);
for k = 1:length(sol)
    u = vals{k};
    if ( ~isreal(u) )
        u = abs(u);
    end
    h(k) = mesh(x{k}, y{k}, u, varargin{:});
    hold on
end
view(0, 90)

% Fix limits (we don't want to do axis equal because of z)
xl = xlim; yl = ylim;
if ( diff(xl) < diff(yl) )
    xlim(xl+(diff(yl)-diff(xl))*[-.5 .5])
else
    ylim(yl+(diff(xl)-diff(yl))*[-.5 .5])
end

if ( ~holdState )
    hold off
end

if ( nargout > 0 )
    varargout = {h};
end

end
