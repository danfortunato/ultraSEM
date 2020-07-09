function c = merge(a, b)
%MERGE   Merge two patch objects.
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
[i1, i2, s1, s2, l2g1, l2g2, scl1, scl2, scl, dom, edges] = intersect(a, b);

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
% - The Dirichlet-to-Neumann maps have their Jacobians factored out.
%   Therefore, we need to multiply the continuity conditions by the
%   multiplication matrices SCL1 and SCL2. The coordinate maps are such
%   that multiplying D2NA/B by SCL1/2 cancels out, so D2NA/B should only be
%   multiplied by SCL1/2, respectively.
A = -( scl2*l2g1*D2Na(s1,s1)*l2g1.' + scl1*l2g2*D2Nb(s2,s2)*l2g2.' );
z = [ scl2*l2g1*D2Na(s1,i1), ...
      scl1*l2g2*D2Nb(s2,i2), ...
      scl2*l2g1*D2Na(s1,end) + scl1*l2g2*D2Nb(s2,end) ];
%    |---------------------- rhs --------------------|

% Select tolerance for least squares minimum norm solver: 
tol = min(max(size(A))*eps(norm(A)), 1e-8);

if ( verLessThan('matlab', '9.3') )
% lsqminnorm() was introduced in MATLAB 2017b. Add code here to make
% ultraSEM *mostly* run on MATLAB 2017a and earlier. 
    dA = A; 
    [U, Sigma, V] = svd( A ); 
    idx = find(diag(Sigma)>tol, 1, 'last');
    ii_remove = idx+1:size(A,2); 
    U(:, ii_remove ) = []; 
    Sigma( ii_remove, : ) = [];
    Sigma( :, ii_remove ) = [];
    V(:, ii_remove) = [];
    S = V * ( Sigma \ ( U'*z ) );
else
% On MATLAB 2017b or later, we can use lsqminnorm. 
% Store the decomposition from lsqminnorm() for reuse in updateRHS().
    dA = matlab.internal.decomposition.DenseCOD(A, tol);
    S = solve(dA, z, false);
end

% Compute new D2N maps:
% - Again, make sure to map local DOFs to shared DOFs on the interface.
Z12 = zeros(numel(i1), numel(i2));
%                                 |--- rhs ----|
D2N = [ D2Na(i1,i1),  Z12,         D2Na(i1,end) ;
        Z12.',        D2Nb(i2,i2), D2Nb(i2,end) ] ...
    + [ D2Na(i1,s1)*l2g1.' ; D2Nb(i2,s2)*l2g2.' ] * S;

% Construct the new patch:
c = ultraSEM.Parent(dom, S, D2N, dA, edges, scl, a, b, ...
                    {i1, s1}, {i2, s2}, l2g1, l2g2, scl1, scl2);

end
