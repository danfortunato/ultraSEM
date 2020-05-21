function f = diffx(f, n)
%DIFFX   Differentiate an ULTRASEM.SOL with respect to x.
%   DIFFX(F) returns an ULTRASEM.SOL representing the derivative of the
%   ULTRASEM.SOL F in its first argument.
%
%   DIFFX(F, N) returns an ULTRASEM.SOL representing the N-th derivative of
%   the ULTRASEM.SOL F in its first argument.
%
% See also DIFFY, DIFF.

% Default to first derivative:
if ( nargin == 1 )
    n = 1;
end

f = diff(f, n, 2);

end
