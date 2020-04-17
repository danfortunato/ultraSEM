function [d, ia, ib] = intersect(S, T)
%INTERSECT   Compute the intersecting patches of two ULTRASEM.DOMAINs.

[d, ia, ib] = intersect(S.domain, T.domain, 'rows', 'stable');

end
