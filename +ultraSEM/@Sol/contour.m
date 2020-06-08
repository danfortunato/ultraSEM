function varargout = contour(sol, varargin)
%CONTOUR   Contour plot of an ULTRASEM.SOL.
%   CONTOUR(SOL) is a contour plot of SOL treating the values of SOL as
%   heights above a plane. A contour plot contains the level curves of SOL
%   for some values V. The values V are chosen automatically.
%
%   CONTOUR(SOL, N) draws N contour lines, choosing the levels
%   automatically.
%   
%   CONTOUR(SOL, V) draws a contour line for each level specified in the
%   vector V. Use CONTOUR(SOL, [V V]) to compute a single contour at the
%   level V.
%
%   [C, H] = CONTOUR(...) returns contour matrix C and a handle H to a
%   contour object. These can be used as inputs to CLABEL. The structure of
%   a contour matrix is described in the help for CONTOURC.
%
%   See also PLOT, SURF, PCOLOR.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

holdState = ishold();
N = 10;
levels = [];

% See if an N was given:
if ( nargin > 1 && isnumeric(varargin{1}) )
    v1 = varargin{1};
    if ( isscalar(v1) )
        N = v1;
    else
        levels = v1;
    end
    varargin(1) = [];
end

% Determine some levels:
if ( isempty(levels) )
    minu = minEst(sol);
    maxu = maxEst(sol);
    levels = linspace(minu, maxu, N);
end

% Loop over the patches:
[x, y] = plotpts(sol);
for k = 1:length(sol)
    u = coeffs2plotvals(sol.coeffs{k});
    if ( ~isreal(u) )
        u = abs(u);
    end
    [c{k}, h(k)] = contour(x{k}, y{k}, u, levels, varargin{:});
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

% Return plot handle if appropriate
if ( nargout == 1 )
    varargout = {h};
elseif ( nargout == 2 )
    varargout = {c, h};
end

end
