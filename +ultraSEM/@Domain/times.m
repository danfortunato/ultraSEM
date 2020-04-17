function varargout = times(varargin)
%.*   Scale an ULTRASEM.DOMAIN.
%   C.*T will scale the ULTRASEM.DOMAIN by C. C must be a scalar.
%
% See also MTIMES.

[varargout{1:nargout}] = mtimes(varargin{:});

end
