function varargout = scatter(sol, varargin)
%SCATTER   Scatter plot of an ULTRASEM.SOL.
%   SCATTER(SOL) displays markers at the tensor product Chebyshev grids for
%   the patches of SOL, colored according to the values of SOL.
%
%   H = SCATTER(...) returns handles to the scatter objects created.
%
%   See also PLOT, MESH.

numPatches = size(sol.domain, 1);
holdState = ishold();

vals = coeffs2vals(sol.coeffs);
[x,y] = getGrid(sol);
for k = 1:numPatches
    h(k) = scatter(x{k}, y{k}, vals{k}, 'EdgeAlpha', 0, varargin{:});
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
