function varargout = plot(varargin)
%PLOT   Plot an ULTRASEM.SOL.
%   PLOT(SOL) gives a surface plot of SOL, the same as SURF(SOL).
%
%   See also SURF, MESH, CONTOUR, PCOLOR.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

[varargout{1:nargout}] = surf(varargin{:});

end
