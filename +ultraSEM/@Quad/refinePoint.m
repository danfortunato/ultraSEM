function [Q, idx] = refinePoint(Q, z)
%REFINEPOINT   Refine an ULTRASEM.QUAD around a point.
%
%   See also REFINE, REFINECORNER.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

loc = find(ismember(Q.v, z, 'rows'));
if ( any(loc) )
    [Q, idx] = refineCorner(Q, loc);
else
    idx = {[1, NaN]};
end

end
