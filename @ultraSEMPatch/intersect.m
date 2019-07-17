function [i1, i2, i4a, i4b, newDom] = intersect(a, b)
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
%   the 'domain' of the intersection. For the time-being this is usually a
%   NaN unless A and B are both rectangular patches whose intersection is
%   also a rectangular patch, i.e.,
%       [ A ][ B ] or [ A ], but not [ A ][ B ] or [ A ]   .
%                     [ B ]               [   ]    [ B    ]
%   This is useful as the indicies of such patches are easier to compute
%   and don't require a call to built-in/intersect() method (which can be
%   expensive).

% Copyright 2018 by Nick Hale and Dan Fortunato.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DEVELOPER NOTES:
% * We currently assume that a mesh has no hanging nodes, so that
%   intersections always occur between entire boundaries.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The new domain will be NaN except in the special case when we are
% merging two rectangluar patches horizontally or vertically.

% Remove corners before we test intersection:
axy = a.xy; bxy = b.xy;
for k = 1:length(a.xy)
    axy{k} = axy{k}(2:end-1,:);
end
for k = 1:length(b.xy)
    bxy{k} = bxy{k}(2:end-1,:);
end
% CAT() is 10x faster than CELL2MAT().
axy = cat(1, axy{:});
bxy = cat(1, bxy{:});

% Convert to single precision to avoid issues with rounding error.
shift = 42;
axy = single(axy + shift);
bxy = single(bxy + shift);

% Determine the common indicies:
if ( ~isnumeric(a.domain) || ~isnumeric(b.domain) )

    [~, i4a, i4b] = intersect(axy, bxy, 'rows', 'stable');
    newDom = ultraSEMDomain([a.domain ; b.domain]);

elseif ( numel(a.domain) == 4 && numel(b.domain) == 4 )
    % Merge two rectangles.

    % We try to avoid calling INTERSECT() if possible (for speed).
    if ( sum(a.domain == b.domain([2 1 3 4])) == 3 )
        % Simple horizontal merge:
        x4 = single(max(a.domain(1), b.domain(1))) + shift;
        i4a = find(axy(:,1) == x4);
        i4b = find(bxy(:,1) == x4);
        % Make sure the ordering of the indicies is correct:
        [~,i1] = sort(axy(i4a,2)); i4a = i4a(i1);
        [~,i2] = sort(bxy(i4b,2)); i4b = i4b(i2);
        newDom = [ min(a.domain(1), b.domain(1)), ...
                   max(a.domain(2), b.domain(2)), ...
                   a.domain(3:4)];

    elseif ( sum(a.domain == b.domain([1 2 4 3])) == 3 )
        % Simple vertical merge:
        y4 = single(max(a.domain(3), b.domain(3))) + shift;
        i4a = find(axy(:,2) == y4);
        i4b = find(bxy(:,2) == y4);
        % Make sure the ordering of the indicies is correct:
        [~,i1] = sort(axy(i4a,1)); i4a = i4a(i1);
        [~,i2] = sort(bxy(i4b,1)); i4b = i4b(i2);
        newDom = [ a.domain(1:2), ...
                   min(a.domain(3), b.domain(3)), ...
                   max(a.domain(4), b.domain(4)) ];

    else
        % Resort to calling intersect:
        [~, i4a, i4b] = intersect(axy, bxy, 'rows', 'stable');
        newDom = ultraSEMDomain([a.domain ; b.domain]);
    end

else
    % Merge two general domains.
    [~, i4a, i4b] = intersect(axy, bxy, 'rows', 'stable');
    newDom = ultraSEMDomain([a.domain ; b.domain]);
end

% We forgot about the "corner" indicies. Add them back.
na = size(a.xy{1},1); ia = unique(floor((i4a-1)/(na-2))+1, 'stable').';
nb = size(b.xy{1},1); ib = unique(floor((i4b-1)/(nb-2))+1, 'stable').';

% Determine what regions of the glue are flipped.
flags1 = find(all(diff(reshape(i4a,na-2,numel(ia))) == -1));
flags2 = find(all(diff(reshape(i4b,nb-2,numel(ib))) == -1));
if ( isempty(i4a) ), flags1 = []; end
if ( isempty(i4b) ), flags2 = []; end
if ( any(flags1) || any(flags2) )
    error('should not need to flip!')
end

% Now create the non-flipped glue indices.
i4a = (1:na).' + (ia-1)*na; i4a = i4a(:);
i4b = (1:nb).' + (ib-1)*nb; i4b = i4b(:);

% i1 and i2 are remaining points (i.e., those not in the intersection).
i1 = (1:(size(axy,1)+2*length(a.xy))).'; i1(i4a) = [];
i2 = (1:(size(bxy,1)+2*length(b.xy))).'; i2(i4b) = [];

% % To "flip" the corresponding glue coefficients, multiply them by
% % [1 -1 1 -1 ...]
% flip1 = ones(size(i4a,1),1);
% flip2 = ones(size(i4b,1),1);
% if ~isempty(flags1), flip1((2:2:na)'+(flags1-1)*na) = -1; end
% if ~isempty(flags2), flip2((2:2:nb)'+(flags2-1)*nb) = -1; end


end
