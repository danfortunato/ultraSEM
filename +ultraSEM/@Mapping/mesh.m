function varargout = mesh(T, varargin)
%MESH   Plot an ULTRASEM.MAPPING as a mesh.
%   MESH(T) plots the ULTRASEM.MAPPING T as a mesh.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

holdState = ishold();

hold on
for k = 1:numel(T)
    % Plot domain:
    vertices = T(k).v;
    h(k,1) = fill(vertices([1:end, 1],1), vertices([1:end, 1],2), ...
        0*vertices([1:end, 1],1), 'FaceAlpha', 0, varargin{:});
end

if ( ~holdState )
    hold off
end

if ( nargout > 0 )
    varargout = {h};
end

end
