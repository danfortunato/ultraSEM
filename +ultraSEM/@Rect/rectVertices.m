function v = rectVertices(R)
%RECTVERTICES   Get the bounding box of an ULTRASEM.RECT.
%   V = RECTVERTICES(R) returns the bounding box of the ULTRASEM.RECT R as
%   a 1x4 vector.
%
%   See also QUADVERTICES.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

v = util.quad2rect(R.v);

end
