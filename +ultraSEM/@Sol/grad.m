function [fx, fy] = grad(f)
%GRAD   Gradient of an ULTRASEM.SOL.
%   [FX FY] = GRAD(F) returns the numerical gradient of the ULTRASEM.SOL F,
%   where FX is the derivative of F in the x direction and FY is the
%   derivative of F in the y direction. Both derivatives are returned as
%   ULTRASEM.SOL objects.
%
%   This is shorthand for GRADIENT(F).
%
% See also GRADIENT.

[fx, fy] = gradient(f);

end
