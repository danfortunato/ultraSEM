function T = removePatch(T, k)
%REMOVEPATCH   Remove patch(es) from an ULTRASEM.DOMAIN.
%   T = REMOVEPATCH(T, K) removes the K-th patch(es) from the
%   ULTRASEM.DOMAIN T and updates the first level mergeIdx appropriately
%   (by introducing NaNs). K may be a vector.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

% Sort for convenience:
k = unique(sort(k(:)));
% Remove specified patches:
T.domain(k,:) = [];

if ( size(T.domain, 1) == 0 )
    % We have removed all the patches!
    T = ultraSEM.Domain();
    return
end

% Remove the corresponding indexes from the first level:
idx1 = T.mergeIdx{1};
for j = 1:numel(k)
    idx1(idx1 == k(j)) = NaN; % Replace missing patches with NaNs.
    idx1(idx1 > k(j)) = idx1(idx1 > k(j)) - 1; % Update.
    k = k - 1;
end
T.mergeIdx{1} = idx1;

end
