function [d, ia, ib] = intersect(S, T)
%INTERSECT   Compute the intersecting patches of two ULTRASEM.DOMAINs.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

[d, ia, ib] = intersect(S.domain, T.domain, 'rows', 'stable');

end
