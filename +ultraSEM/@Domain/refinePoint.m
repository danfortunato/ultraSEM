function [T, newIdx] = refinePoint(T, z, m)
%REFINEPOINT   Refine an ULTRASEM.DOMAIN around a point.
%
%   See also REFINE.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

newIdx = {};
if ( nargin < 3 ), m = 1; end
if ( m == 0 ), return, end

if ( numel(T) > 1 )
    for k = numel(T):-1:1
        T(k) = refinePoint(T(k), z, m);
    end
    return 
end
n = numel(T.domain);
newDom = cell(n,1);
idx = cell(n,1);
s = [0, 0, 0];

for k = 1:n
    [newDom{k}, newIdx] = refinePoint(T.domain(k), z);
    
    if ( numel( newDom{k} ) == 1 )
          newIdx = {[1, NaN], [1, NaN], [1, NaN]};  
    end
    
    if ( numel(newIdx) == 2 )
        newIdx = [newIdx, {[1, NaN]}];
    end
    
    % Append new merge indices (with appropriate shifts):
    for j = 1:numel(newIdx)
        idx{k,j} = newIdx{j} + s(j);
        s(j) = max(idx{k,j}(:));
    end
    
end

% TODO: Lazy. Need to do this for general number of levels.
if ( size(idx, 2) == 1 )
    idx = {};
elseif ( size(idx, 2) == 2 )
    idx = {vertcat(idx{:,1}), vertcat(idx{:,2})};
else
    idx = {vertcat(idx{:,1}), vertcat(idx{:,2}), vertcat(idx{:,3})};
end

for k = numel(idx):-1:1
    if (all(isnan(idx{k}(:,2))))
        idx(k) = [];
    end
end

T.domain = vertcat(newDom{:});
T.mergeIdx = [idx, T.mergeIdx];

[T, newIdx] = refinePoint(T, z, m-1);

% Stacked
% newIdx = {};
% if ( nargin < 3 ), m = 1; end
% if ( m == 0 ), return, end
% 
% if ( numel(T) > 1 )
%     for k = 1:numel(T)
%         T(k) = refinePoint(T(k), z, m);
%     end
%     return 
% end
% 
% prevDom = T.domain;
% if ( ~isa(T.domain(1), 'ultraSEM.Domain') )
%     dom = {};
%     for k = 1:numel(T.domain)
%         dom{k} = ultraSEM.Domain(T.domain(k));
%     end
%     T.domain = vertcat(dom{:});
% end
% 
% newDom = {};
% for k = 1:numel(T.domain)
%     [T.domain(k), newIdx] = refinePoint(prevDom(k), z);
%     T.domain(k).mergeIdx = newIdx;
% end
% 
% [T, newIdx] = refinePoint(T, z, m-1);

end