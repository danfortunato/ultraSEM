function c = centroid(Q)
%CENTROID   Compute the centroid of an ULTRASEM.MAPPING.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

[xmid, ymid] = centroid(polyshape(Q.v, 'Simplify', false));
c = [xmid ; ymid];

end
