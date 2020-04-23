function [i1, i2, i4a, i4b, l2g1, l2g2, dom, edgesAB] = intersect(a, b)
%INTERSECT   Compute the indices of the glue between two patches.
%   [I1, I2, I4A, I4B, L2G1, L2G2] = INTERSECT(A, B) returns the indices
%   of the glue w.r.t. A.edges and B.edges of two patches A and B. For
%   consistency with the paper by Martinsson:
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
%   the shared interface.
%
%   [I1, I2, I4A, I4B, L2G1, L2G2, DOM, EDGESAB] = INTERSECT(A, B) returns
%   also the domain and the edges of the intersection.

% Copyright 2018 by Nick Hale and Dan Fortunato.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DEVELOPER NOTES:
% * We currently assume that a mesh has no hanging nodes, so that
%   intersections always occur between entire boundaries.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

% Build interpolation/truncation operators for p-adaptivity.
if ( isempty(i4a) )
    % Force 0x1 rather than empty. (Not sure why this is required.)
    int1 = ones(0,1);
    int2 = ones(0,1);
else
    % Get the function to determine the polynomial degree on the interface.
    pref = ultraSEM.Pref();

    blocks1 = cell(numel(iA),1);
    blocks2 = cell(numel(iB),1);
    for k = 1:numel(iA)
        % Set block (k,k) of int1 and int2
        pmax = pref.interfaceDegree(pA(iA(k)), pB(iB(k)));
        blocks1{k} = speye(pmax, pA(iA(k)));
        blocks2{k} = speye(pmax, pB(iB(k)));
    end
    int1 = blkdiag(blocks1{:});
    int2 = blkdiag(blocks2{:});
end

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