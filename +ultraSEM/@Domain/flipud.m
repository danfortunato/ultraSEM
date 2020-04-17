function T = flipud(T)
%FLIPUD   Flip an ULTRASEM.DOMAIN vertically about its center.

miny = min(T.domain(:,3));
maxy = max(T.domain(:,4));
T.domain(:,[3 4]) = miny + maxy - T.domain(:,[4 3]);

end
