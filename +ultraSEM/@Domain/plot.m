function varargout = plot(T, varargin)
%PLOT   Plot an ULTRASEM.DOMAIN.
%   PLOT(T) plots the ULTRASEM.DOMAIN T and indicates the patch numbering
%   for the first layer.
%
%   PLOT(T, C), where C is a single character string chosen from the list
%   'r', 'g', 'b', 'c', 'm', 'y', 'w', 'k', or an RGB row vector triple,
%   [r g b], fills the domain with the constant specified color.
%
%   PLOT(T, C, 'PropertyName', PropertyValue, ...) allows adjusting of the
%   plot in any way allowed by the built-in FILL() method.
%
%   H = PLOT(T, ...) returns a figure handle of the form returned by
%   H = FILL(...), where FILL() is the built-in MATLAB method.

if ( numel(T) > 1 ) 
    holdState = ishold();
    for k = 1:numel(T)
        plot(T(k), varargin{:}); hold on
    end
    if ( ~holdState )
        hold off
    end
else
    [varargout{1:nargout}] = plot(T.domain, varargin{:});
end

end
