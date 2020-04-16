function varargout = pcolor(sol, varargin)
%PCOLOR   Pseudocolor (checkerboard) plot of an ULTRASEM.SOL.
%   PCOLOR(SOL) is a pseudocolor or "checkerboard" plot of the solution
%   SOL.
%
%   H = PCOLOR(...) returns a handle to a SURFACE object.
%
% See also PLOT, SURF, CONTOUR.

holdState = ishold();

[x,y] = plotpts(sol);

for k = 1:length(sol)
    u = coeffs2plotvals(sol.u{k});
    if ( ~isreal(u) )
        u = abs(u);
    end
    h(k) = surface(x{k}, y{k}, 0*u, u, 'EdgeAlpha', 1, varargin{:});
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
