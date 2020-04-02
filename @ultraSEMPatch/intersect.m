function [i1, i2, i4a, i4b, flip1, flip2, newDom, edgesAB] = intersect(a, b)
%INTERSECT   Compute the indicies of the glue between two patches.
%   [I1, I2, I4A, I4B, FLIP1, FLIP2] = INTERSECT(A, B) returns the indicies
%   of the glue w.r.t. A.xy and B.xy of two patches A and B. For
%   consistency with the paper by Martinsson:
%       I1 : Indicies of points A.xy which are not in B.xy
%       I2 : Indicies of points B.xy which are not in A.xy
%       I4A: Indicies of A.xy which are in the intersection
%       I4B: Indicies of B.xy which are in the intersection
%
%   If the boundaries being merged contain flipped regions, so that the
%   glue indicies are ordered differently local to each patch, then FLIP1
%   and FLIP2 encode how to "flip" the shared coefficients for each patch.
%
%   [I1, I2, I4A, I4B, FLIP1, FLIP2, DOM] = INTERSECT(A, B) returns also
%   the 'domain' of the intersection. 

% Copyright 2018 by Nick Hale and Dan Fortunato.

% TODO: Update help text.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DEVELOPER NOTES:
% * We currently assume that a mesh has no hanging nodes, so that
%   intersections always occur between entire boundaries.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The new domain will be NaN except in the special case when we are
% merging two rectangluar patches horizontally or vertically.
newDom = NaN;

% Above to eventually be replaced by:
edgesA = a.edges;
edgesB = b.edges;

% Check for intersecting edges (with a tolerance):
[iA1,iB1] = intersectTol(edgesA(:,[1 2 3 4]), edgesB(:,1:4), 1e-12);
% Must check in both directions:
[iA2,iB2] = intersectTol(edgesA(:,[3 4 1 2]), edgesB(:,1:4), 1e-12);
iA = [iA1 ; iA2];
iB = [iB1 ; iB2];

% Determine indices corresponding to intersecting edges:
pA = edgesA(:,5);
ppA = cumsum([0 ; pA]);
i4a = (ppA(iA)+(1:pA(iA))).'; i4a = i4a(:);
pB = edgesB(:,5);
ppB = cumsum([0 ; pB]);
i4b = (ppB(iB)+(1:pB(iB))).'; i4b = i4b(:);

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

edgesA(iA,:) = [];
edgesB(iB,:) = [];
edgesAB = [edgesA ; edgesB];

try
    newDom = ultraSEMDomain([a.domain, b.domain]); % Only for debugging. Remove.
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