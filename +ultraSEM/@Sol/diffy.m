function f = diffy(f)
%DIFFY   Differentiate an ULTRASEM.SOL with respect to y.
%   DIFFX(F) returns an ULTRASEM.SOL representing the derivative of the
%   ULTRASEM.SOL F in its second argument.
%
%   DIFFX(F, N) returns an ULTRASEM.SOL representing the N-th derivative of
%   the ULTRASEM.SOL F in its second argument.
%
%   See also DIFFX, DIFF.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

% Default to first derivative:
if ( nargin == 1 )
    n = 1;
end

f = diff(f, n, 1);

end
