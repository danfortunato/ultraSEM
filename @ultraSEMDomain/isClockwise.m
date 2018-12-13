function out = isClockwise( v )
%ISCLOCKWISE   Determines whether a polgon is clokwise orientated.
%   ISCLOCKWISE(V) returns true if the polygon defined by the vertices V is
%   clockwise oriented, or not. It does this by computing the signed area
%   of the polygon. If the polygon is degenerate (i.e., the area is zero)
%   then NaN is returned.

% Given a polygon with vertices (x(1),y(1)),...,(x(n),y(n)), the signed
% area is 1/2 (sum (x(i)*y(i+1) - x(i+1)*y(i) )
v1 = v(:,1); v2 = v(:,2);
signedArea = sum( v1.*circshift(v2,1) - circshift(v1,1).*v2 ) / 2;

% If the signed area is positive then the polygon has a clockwise
% orientation:
out = signedArea > 0;

% Return NaN if we encounter a degenerate polygon.
if ( signedArea == 0 )
    out = NaN;
end

end
