function c = merge(a, b)
%MERGE   Merge two ultraSEMPatch objects.
%   C = MERGE(A, B) returns a patch C formed by merging the two patches A
%   and B. Typically A and B will be adjacent and any common edges will be
%   eliminated by enforcing continuity and continuity of the derivative
%   across the boundary.

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
[i1, i2, s1, s2, l2g1, l2g2, dom, edges] = intersect(a, b);

% Extract D2N maps:
D2Na = a.D2N; D2Nb = b.D2N;
% and discard from children
% a.D2N = []; b.D2N = []; % Keep so that we can update the RHS efficiently(?)

% Compute new solution operator:
% - In the case of cross points, the operator we have to invert is
%   rank deficient by one. The least squares solution is correct, though.
% - Make sure to map local DOFs on A and B to global DOFs on the shared
%   interface. The L2G maps are responsible for any transformations needed
%   to achieve this (e.g. flipping, interpolation, etc.).
S = lsqminnorm( -( l2g1*D2Na(s1,s1)*l2g1.' + l2g2*D2Nb(s2,s2)*l2g2.' ), ...
                 [ l2g1*D2Na(s1,i1), ...
                   l2g2*D2Nb(s2,i2), ...
                   l2g1*D2Na(s1,end) + l2g2*D2Nb(s2,end) ]);
%                 |----------------- rhs ---------------|

% Compute new D2N maps:
% - Again, make sure to map local DOFs to shared DOFs on the interface.
Z12 = zeros(numel(i1), numel(i2));
%                                 |--- rhs ----|
D2N = [ D2Na(i1,i1),  Z12,         D2Na(i1,end) ;
        Z12.',        D2Nb(i2,i2), D2Nb(i2,end) ] ...
    + [ D2Na(i1,s1)*l2g1.' ; D2Nb(i2,s2)*l2g2.' ] * S;

% Construct the new patch:
c = ultraSEMParent(dom, S, D2N, edges, a, b, {i1, s1}, {i2, s2}, l2g1, l2g2);
end
