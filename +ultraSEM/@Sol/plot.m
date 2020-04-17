function varargout = plot(varargin)
%PLOT   Plot an ULTRASEM.SOL.
%   PLOT(SOL) gives a surface plot of SOL, the same as SURF(SOL).
%
% See also SURF, MESH, CONTOUR, PCOLOR.

[varargout{1:nargout}] = surf(varargin{:});

end
