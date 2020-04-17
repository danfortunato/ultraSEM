function v = rectVertices(R)
%RECTVERTICES   Get the bounding box of an ULTRASEM.RECT.
%   V = RECTVERTICES(R) returns the bounding box of the ULTRASEM.RECT R as
%   a 1x4 vector.
%
% See also QUADVERTICES.

v = util.quad2rect(R.v);

end
