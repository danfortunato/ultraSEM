function varargout = surf(sol, varargin)
%SURF   3-D colored surface of an ULTRASEM.SOL.
%   SURF(SOL) plots a colored parametric surface whose height and color is
%   defined by the values of SOL.
%
%   SURF(..., 'PropertyName', PropertyValue, ...) sets the value of the
%   specified surface property. Multiple property values can be set with a
%   single statement.
%
%   H = SURF(...) returns a handle to a surface plot object.
%
% See also PLOT, CONTOUR, MESH, PCOLOR.

holdState = ishold();

[x, y] = plotpts(sol);
if ( ~iscell(x) )
    x = {x};
    y = {y};
end

for k = 1:length(sol)
    u = coeffs2plotvals(sol.u{k});
    if ( ~isreal(u) )
        u = abs(u);
    end
    h(k) = surf(x{k}, y{k}, u, 'EdgeAlpha', 1, varargin{:});
    hold on
end
shading interp
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
