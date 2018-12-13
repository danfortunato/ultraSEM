function [i1, i2, i4a, i4b, newDom] = intersect(a, b)
%INTERSECT   Compute the indicies of the glue between two patches.
%   [I1, I2, I4A, I4B] = INTERSECT(A, B) returns the indicies of the glue
%   w.r.t. A.xy and B.xy of two patches A and B. For consistency with the
%   paper by Martinsson:
%       I1 : Indicies of points A.xy which are not in B.xy
%       I2 : Indicies of points B.xy which are not in A.xy
%       I4A: Indicies of A.xy which are in the intersection
%       I4B: Indicies of B.xy which are in the intersection
%
%   [I1, I2, I4A, I4B, DOM] = INTERSECT(A, B) returns also the 'domain' of the
%   intersection. For the time-being this is usually a NaN unless A and B are
%   both rectangular patches whose intersection is also a rectangular patch,
%   i.e.,
%       [ A ][ B ] or [ A ], but not [ A ][ B ] or [ A ]   .
%                     [ B ]               [   ]    [ B    ]
%   This is useful as the indicies of such patches are easier to compute and
%   don't require a call to built-in/intersect() method (which can be
%   expensive).

% Copyright 2018 by Nick Hale and Dan Fortunato.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DEVELOPER NOTES:
% * We currently assume that a mesh has no hanging nodes, so that
%   intersections always occur between entire boundaries.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The new domain will be NaN except in the special case when we are
% merging two rectangluar patches horizontally or vertically.
newDom = NaN;

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

elseif ( numel(a.domain) == 4 && numel(b.domain) == 4 )
    % Merge two rectangles.

    % We try to avoid calling INTERSECT() if possible (for speed).
    if ( sum(a.domain == b.domain([2 1 3 4])) == 3 )
        % Simple horizontal merge:
        x4 = single(max(a.domain(1), b.domain(1))) + shift;
        i4a = find(axy(:,1) == x4);
        i4b = find(bxy(:,1) == x4);
        newDom = [ min(a.domain(1), b.domain(1)), ...
                   max(a.domain(2), b.domain(2)), ...
                   a.domain(3:4)];

    elseif ( sum(a.domain == b.domain([1 2 4 3])) == 3 )
        % Simple vertical merge:
        y4 = single(max(a.domain(3), b.domain(3))) + shift;
        i4a = find(axy(:,2) == y4);
        i4b = find(bxy(:,2) == y4);
        newDom = [ a.domain(1:2), ...
                   min(a.domain(3), b.domain(3)), ...
                   max(a.domain(4), b.domain(4)) ];

    else
        % Resort to calling intersect:
        [~, i4a, i4b] = intersect(axy, bxy, 'rows', 'stable');

    end

else
    % Merge two general domains.
    [~, i4a, i4b] = intersect(axy, bxy, 'rows', 'stable');
end

% i1 and i2 are remaining points (i.e., those not in the intersection.
i1 = (1:size(axy,1)).'; i1(i4a) = [];
i2 = (1:size(bxy,1)).'; i2(i4b) = [];

end
