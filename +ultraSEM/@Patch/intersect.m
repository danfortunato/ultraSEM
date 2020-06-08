function [i1, i2, i4a, i4b, l2g1, l2g2, scl1, scl2, sclAB, dom, edgesAB] = intersect(a, b)
%INTERSECT   Compute the indices of the glue between two patches.
%   [I1, I2, I4A, I4B, L2G1, L2G2, SCL1, SCL2, SCLAB] = INTERSECT(A, B)
%   returns the indices of the glue w.r.t. A.edges and B.edges of two
%   patches A and B. For consistency with the paper by Martinsson:
%
%       I1 : Indices of DOFs on A.edges which are not on B.edges
%       I2 : Indices of DOFs on B.edges which are not on A.edges
%       I4A: Indices of DOFs on A.edges which are in the intersection
%       I4B: Indices of DOFs on B.edges which are in the intersection
%
%   If the boundaries being merged contain flipped regions (so that DOFs
%   are ordered differently local to each patch) or regions with different
%   polynomial degrees on either side of an interface (so that DOFs must be
%   interpolated or restricted to the same space), then L2G1 and L2G2
%   encode how to map from local DOFs in A and B to global DOFs along
%   the shared interface. The matrices SCL1 and SCL2 are multiplication
%   matrices for the algebraic expressions (e.g., Jacobians) that have been
%   factored out of the Dirichlet-to-Neumann maps for patches A and B.
%   SCLAB is a cell array of scalars and/or function handles defining those
%   algebraic expressions for the parent's edges.
%
%   [I1, I2, I4A, I4B, L2G1, L2G2, SCL1, SCL2, SCLAB, DOM, EDGESAB] = INTERSECT(A, B)
%   returns also the domain and the edges of the intersection.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DEVELOPER NOTES:
% * We currently assume that a mesh has no hanging nodes, so that
%   intersections always occur between entire boundaries.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pref = ultraSEM.Pref();

% The new domain will be NaN except in the special case when we are
% merging two rectangluar patches horizontally or vertically.
dom = NaN;

% Above to eventually be replaced by:
edgesA = a.edges;
edgesB = b.edges;

% Check for intersecting edges (with a tolerance):
[iA1,iB1] = intersectTol(edgesA(:,[1 2 3 4]), edgesB(:,1:4), 1e-12);
% Must check in both directions:
[iA2,iB2] = intersectTol(edgesA(:,[3 4 1 2]), edgesB(:,1:4), 1e-12);
iA = [iA1 ; iA2];
iB = [iB1 ; iB2];
assert(numel(iA) == numel(iB), 'Intersection failed.');

% Determine indices corresponding to intersecting edges:
pA = edgesA(:,5);
ppA = cumsum([0 ; pA]);
i4a = [];
for k = 1:numel(iA)
    i4a = [i4a ; ppA(iA(k)) + (1:pA(iA(k))).'];
end
pB = edgesB(:,5);
ppB = cumsum([0 ; pB]);
i4b = [];
for k = 1:numel(iB)
    i4b = [i4b ; ppB(iB(k)) + (1:pB(iB(k))).'];
end

% i1 and i2 are remaining points (i.e., those not in the intersection).
i1 = (1:sum(pA)).'; i1(i4a) = [];
i2 = (1:sum(pB)).'; i2(i4b) = [];

% TODO: Deal with flipping?
flip1 = ones(size(i4a));
flip2 = ones(sum(pB(iB1)),1);
for k = 1:numel(iB2)
    e = ones(pB(iB2(k)),1); e(2:2:end) = -1;
    flip2 = [flip2 ; e];
end
if ( isempty(flip1) )
    % Force 0x1 rather than empty. (Not sure why this is required.)
    flip1 = ones(0,1);
    flip2 = ones(0,1);
end

%% Construct operators for p-adaptivity and Jacobian scaling.

sclA = a.D2N_scl;
sclB = b.D2N_scl;

intBlocks1 = cell(numel(iA),1); intBlocks2 = cell(numel(iB),1);
sclBlocks1 = cell(numel(iA),1); sclBlocks2 = cell(numel(iB),1);

for k = 1:numel(iA)
    % Determine the polynomial degree on the interface.
    n = pref.interfaceDegree(pA(iA(k)), pB(iB(k)));

    %% Build interpolation/truncation operators for p-adaptivity.
    % Block (k,k) interpolates from patch A/B to interface iA(k)/iB(k).
    intBlocks1{k} = speye(n, pA(iA(k)));
    intBlocks2{k} = speye(n, pB(iB(k)));

    %% Compute the Jacobian scaling for patch A at interface iA(k).
    s = sclA{iA(k)};
    if ( isnumeric(s) && isscalar(s) )
        sclBlocks1{k} = s * speye(n);
    else
        % We need to build a multiplication matrix for the scaling, so
        % evaluate at Chebyshev points to convert to coefficients.
        x = util.chebpts(n);
        cfsA = util.vals2coeffs(s(x));
        sclBlocks1{k} = util.multmat(n, cfsA);
    end

    %% Compute the Jacobian scaling for patch B at interface iB(k).
    s = sclB{iB(k)};
    if ( isnumeric(s) && isscalar(s) )
        sclBlocks2{k} = s * speye(n);
    else
        % We need to build a multiplication matrix for the scaling, so
        % evaluate at Chebyshev points to convert to coefficients.
        x = util.chebpts(n);
        cfsB = util.vals2coeffs(s(x));
        sclBlocks2{k} = util.multmat(n, cfsB);
    end
end

% Assemble the operators as block-diagonal matrices.
int1 = blkdiag(intBlocks1{:}); int2 = blkdiag(intBlocks2{:});
scl1 = blkdiag(sclBlocks1{:}); scl2 = blkdiag(sclBlocks2{:});

% Concatenate the scaling functions for the parent's edges.
if ( ~isempty(iA) )
    sclA(iA) = [];
    sclB(iB) = [];
end
sclAB = [sclA ; sclB];

% Local-to-global maps encode interpolation and flipping.
l2g1 = int1 .* flip1.';
l2g2 = int2 .* flip2.';

edgesA(iA,:) = [];
edgesB(iB,:) = [];
edgesAB = [edgesA ; edgesB];

try
    dom = ultraSEM.domain([a.domain, b.domain]); % Only for debugging. Remove.
end

end


function [AI, BI] = intersectTol(A, B, tol)
% Borrowed from https://www.mathworks.com/matlabcentral/...
% answers/444501-how-to-use-intersect-command-with-a-tolerance-value
   n = size(A,1);
   M  = zeros(n,1);
   % Collect the index of the first occurrence in B for every A:
   for k = 1:n
      dist = sum(abs(A(k,:) - B), 2);    % 1-norm
      idx  = find(dist < tol, 1);        % Absolute tolerance
      % iddx = find(dist ./ A(iA) < tol, 1);  % Relative tolerance
      if ( ~isempty(idx) )
         M(k) = idx;
      end
   end
   AI = find(M);
   BI = M(AI); 
end