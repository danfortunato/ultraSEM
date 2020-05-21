function f = diffx(f)
%DIFFX   Differentiate an ULTRASEM.SOL with respect to x.
%   DIFFX(F) returns an ULTRASEM.SOL representing the derivative of the
%   ULTRASEM.SOL F in its first argument.
%
% See also DIFFY.

f = diff(f, 2);

end
