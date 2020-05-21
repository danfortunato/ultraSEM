function f = diffy(f)
%DIFFY   Differentiate an ULTRASEM.SOL with respect to y.
%   DIFFX(F) returns an ULTRASEM.SOL representing the derivative of the
%   ULTRASEM.SOL F in its second argument.
%
% See also DIFFX.

f = diff(f, 1);

end
