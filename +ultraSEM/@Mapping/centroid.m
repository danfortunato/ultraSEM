function c = centroid(Q)
%CENTROID   Compute the centroid of an ULTRASEM.MAPPING.

[xmid, ymid] = centroid(polyshape(Q.v, 'Simplify', false));
c = [xmid ; ymid];

end
