function varargout = times(varargin)
%.*   Scale an ULTRASEM.SOL.
%   C.*T will scale the ULTRASEM.SOL by C. C must be a scalar.
%
%   See also MTIMES.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

[varargout{1:nargout}] = mtimes(varargin{:});

end
