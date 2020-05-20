function varargout = times(varargin)
%.*   Scale an ULTRASEM.SOL.
%   C.*T will scale the ULTRASEM.SOL by C. C must be a scalar.
%
% See also MTIMES.

[varargout{1:nargout}] = mtimes(varargin{:});

end
