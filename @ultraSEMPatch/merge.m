function c = merge(a, b)
%MERGE   Merge two ultraSEMPatch objects.
%   C = MERGE(A, B) returns a patch C formed by merging the two patches A
%   and B. Typically A and B will be adjacent and any common grid points
%   (i.e., those in A.xy and B.xy) will be eliminated by enforcing
%   continuity and continuity of the derivative across the boundary.

% Copyright 2018 by Nick Hale and Dan Fortunato.

% Parse inputs:
if ( nargin == 0 )
    c = [];
    return
elseif ( nargin == 1 )
    c = a;
    return
end

% Compute the indices of intersecting points in a and b.
[i1, i2, s1, s2, flip1, flip2, dom] = intersect(a, b);

% Extract D2N maps:
D2Na = a.D2N; D2Nb = b.D2N;
% and discard from children
% a.D2N = []; b.D2N = []; % Keep so that we can update the RHS efficiently(?)

% Compute new solution operator:
% - In the case of cross points, the operator we have to invert is
% rank deficient by one. The least squares solution is correct, though.
% - Make sure to "flip" the shared DOFs on the interface.
S = lsqminnorm( -( flip1.*D2Na(s1,s1).*flip1.' + flip2.*D2Nb(s2,s2).*flip2.' ), ...
                 [ flip1.*D2Na(s1,i1), ...
                   flip2.*D2Nb(s2,i2), ...
                   flip1.*D2Na(s1,end) + flip2.*D2Nb(s2,end) ]);
%                 |----------------- rhs -------------------|

% Compute new D2N maps:
% - Again, make sure to "flip" the shared DOFs on the interface.
Z12 = zeros(numel(i1), numel(i2));
%                                 |--- rhs ----|
D2N = [ D2Na(i1,i1),  Z12,         D2Na(i1,end) ;
        Z12.',        D2Nb(i2,i2), D2Nb(i2,end) ] ...
    + [ D2Na(i1,s1).*flip1.' ; D2Nb(i2,s2).*flip2.' ] * S;

% Construct the new patch:
xy = mergeBdy(a.xy, b.xy, i1, i2);
c = ultraSEMParent(dom, S, D2N, xy, a, b, {i1, s1}, {i2, s2}, flip1, flip2);

end

function xy = mergeBdy(axy, bxy, i1, i2)
%MERGBDY   Merge the boundary nodes of two patches.
%   XY = MERGEBDY(AXY, BXY, I1, I2) merges boundaries containing nodes
%   specified by I1 and I2. We currently assume that the mesh has no
%   hanging nodes, so that entire boundaries are always merged with entire
%   boundaries. If for some reason I1 or I2 contain only subset of nodes
%   for a given boundary, we merge all the nodes for that boundary.

% Map global indicies to indicies specifying which boundaries to merge:
na = size(axy{1},1); ia = unique(floor((i1-1)/na)+1, 'stable');
nb = size(bxy{1},1); ib = unique(floor((i2-1)/nb)+1, 'stable');
xy = [axy(ia) ; bxy(ib)];

end