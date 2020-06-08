function T = fliplr(T)
%FLIPLR   Flip an ULTRASEM.DOMAIN horizontally about its center.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

minx = min(T.domain(:,1));
maxx = max(T.domain(:,2));
T.domain(:,[1 2]) = minx + maxx - T.domain(:,[2 1]);

end
