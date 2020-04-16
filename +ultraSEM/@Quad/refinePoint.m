function [Q, idx] = refinePoint(Q, z)
%REFINEPOINT   Refine an ULTRASEM.QUAD around a point.

loc = find(ismember(Q.v, z, 'rows'));
if ( any(loc) )
    [Q, idx] = refineCorner(Q, loc);
else
    idx = {[1, NaN]};
end

end
