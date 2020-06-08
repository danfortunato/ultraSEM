function T = transpose(T)
%.'   Rotate an ULTRASEM.DOMAIN counterclockwise by 90 degrees.
%
%   See also CTRANSPOSE, ROT90.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

T = rot90(T, -1);

end
