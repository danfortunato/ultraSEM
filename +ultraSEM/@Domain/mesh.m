function varargout = mesh(T, varargin)
%MESH   Plot an ULTRASEM.DOMAIN as a mesh.
%   MESH(T) plots the ULTRASEM.DOMAIN T.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

if ( numel(T) > 1 ) 
    holdState = ishold();
    for k = 1:numel(T)
        mesh(T(k), varargin{:}); hold on
    end
    if ( ~holdState )
        hold off
    end
else
    [varargout{1:nargout}] = mesh(T.domain, varargin{:});
end

end
