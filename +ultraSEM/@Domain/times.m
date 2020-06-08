function varargout = times(varargin)
%.*   Scale an ULTRASEM.DOMAIN.
%   C.*T will scale the ULTRASEM.DOMAIN by C. C must be a scalar.
%
%   See also MTIMES.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

[varargout{1:nargout}] = mtimes(varargin{:});

end
