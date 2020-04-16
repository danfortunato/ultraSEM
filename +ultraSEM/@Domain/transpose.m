function T = transpose(T)
%.'   Rotate an ULTRASEM.DOMAIN counterclockwise by 90 degrees.
%
% See also ROT90.

T = rot90(T, -1);

end
