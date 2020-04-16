function [T, newIdx] = refinePoint(T, z, m)
%REFINEPOINT   Refine an ULTRASEM.DOMAIN around a point.

newIdx = {};
if ( nargin < 3 ), m = 1; end
if ( m == 0 ), return, end

if ( numel(T) > 1 )
    for k = 1:numel(T)
        T(k) = refinePoint(T(k), z, m);
    end
    return 
end
n = numel(T.domain);
newDom = cell(n,1);
idx = cell(n,1);
s = [0, 0];
for k = 1:n
    [newDom{k}, newIdx] = refinePoint(T.domain(k), z);
    if ( numel( newIdx ) == 1 )
        newIdx = {[1, NaN], [1, NaN]};
    end
    idx{k,1} = newIdx{1} + s(1);
    idx{k,2} = newIdx{2} + s(2);
    s = [max(idx{k,1}(:)), max(idx{k,2}(:))];
end
idx = {vertcat(idx{:,1}), vertcat(idx{:,2})};

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
