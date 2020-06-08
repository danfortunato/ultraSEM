function varargout = plot(T, varargin)
%PLOT   Plot the mapped domain and grid.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

% Choose a color:
if ( nargin > 1 && ischar(varargin{1}) && ...
        ~isempty(regexp( varargin{1}, '[bgrcmykw]', 'match')) )
    col = varargin{1};
    varargin(1) = [];
elseif ( nargin > 1 && isnumeric(varargin{1}) )
    col = varargin{1};
    varargin(1) = [];
else
    col = rand(1, 3);
end

alpha = 0.25;
if ( size(col,2) == 4 )
    alpha = col(4);
    col = col(1:3);
end
plotPts = false;
n = 21;
holdState = ishold();

if ( nargin > 1 && ~isempty(varargin) && ...
        ~isempty(regexp(varargin{1}, '[.ox+*sdv^<>ph]', 'match')) && ...
        ~any(strcmpi(varargin{1}, {'LineStyle', 'LineWidth'})) )
    plotPts = true;
    marker = varargin{1};
    varargin(1) = [];
end

for k = 1:numel(T)
    % Plot domain:
    vertices = T(k).v;
    h(k,1) = fill(vertices([1:end, 1],1), vertices([1:end, 1],2), col, ...
        'FaceAlpha', alpha, varargin{:}); hold on %#ok<AGROW>

    if ( plotPts )
        % Plot the grid:
        [xx, yy] = util.chebpts2( n, n, T(k) );
        plot(xx, yy, marker, varargin{:})
    end

    % Add text to center of patch:
    if ( numel(T) < 100 ) % Don't add text on large meshes
        c = centroid(T(k));
        text(c(1), c(2), int2str(k), 'HorizontalAlignment', 'center')
    end

end

if ( ~holdState )
    hold off
end

if ( nargout > 0 )
    varargout = {h};
end

end
