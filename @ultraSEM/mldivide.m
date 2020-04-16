function varargout = mldivide(varargin)
%MLDIVIDE   Solve an ULTRASEM system.
%   MLDIVIDE is a shorthand for SOLVE().
%
% See also SOLVE.

% MLDIVIDE() is simply a wrapper for SOLVE():
[varargout{1:nargout}] = solve(varargin{:});

end
