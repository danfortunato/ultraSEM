function [fx, fy] = gradient(f)
%GRADIENT   Gradient of an ULTRASEM.SOL.
%   [FX FY] = GRADIENT(F) returns the numerical gradient of the
%   ULTRASEM.SOL F, where FX is the derivative of F in the x direction and
%   FY is the derivative of F in the y direction. Both derivatives are
%   returned as ULTRASEM.SOL objects.
%
% See also GRAD.

fx = diff(f, 1, 2);
fy = diff(f, 1, 1);

end
