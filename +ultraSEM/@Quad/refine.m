function [Q, mergeIdx] = refine(Q, m)
%REFINE   Refine an ULTRASEM.QUAD.
%
%   See also REFINECORNER, REFINEPOINT.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

% Parse inputs:
mergeIdx = {}; 
if ( nargin == 1 ), m = 1; end
if ( m == 0 ), return, end

n = numel(Q); % Number of Quads

% Initialize:
v = cell(4*n, 1);
idx1 = cell(n,1);
idx2 = cell(n,1);

% Divide into four using centroid and centre of edges:
for k = 1:n
    vk = Q(k).v;
    c = centroid(Q(k)).';
    vnew = (vk([1 2 3 4],:) + vk([2 3 4 1],:))/2;
    v(4*(k-1)+(1:4)) = ...
         {[vk(1,:) ; vnew(1,:) ; c ; vnew(4,:)];
          [vnew(1,:) ; vk(2,:) ; vnew(2,:) ; c];
          [c ; vnew(2,:) ; vk(3,:) ; vnew(3,:)];
          [vnew(4,:) ; c ; vnew(3,:) ; vk(4,:)]};
    idx1{k} = [1 2 ; 3 4] + 4*(k-1);
    idx2{k} = [1 2] + 2*(k-1);
end

% Assemble Quad:
mergeIdx = {vertcat(idx1{:}), vertcat(idx2{:})};
Q = ultraSEM.Quad(v);           

% Recurse for more refinement:
[Q, idxNew] = refine(Q, m-1);

% Append merge info:
mergeIdx = [idxNew, mergeIdx];

end
