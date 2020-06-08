function H = oldMerge(F, G, varargin)
%OLDMERGE   Merge two or more ULTRASEM.DOMAINs (old-style merge).
%
%   See also MERGE.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

if ( nargin == 1 )
    H = F;
    return
end

% Merge multiple pieces:
if ( nargin > 2 )
    H = merge(merge(F,G), varargin{:});
    return
end
% Nothing to do here:
if ( isempty(F) )
    H = G;
    return
elseif ( isempty(G) )
    H = F;
    return
end

treeSizeF = size(F.mergeIdx, 2);
treeSizeG = size(G.mergeIdx, 2);
% The output will have one more level than F or G.
treeSizeH = max(treeSizeF, treeSizeG) + 1;

% Pad F or G so that they have the same number of levels:
if ( treeSizeF < treeSizeG )
    F.mergeIdx(end+1:treeSizeG) = ...
        repmat({[1, NaN]}, 1, treeSizeG-treeSizeF);
elseif ( treeSizeF > treeSizeG )
    G.mergeIdx(end+1:treeSizeF) = ...
        repmat({[1, NaN]}, 1, treeSizeF-treeSizeG);
end

% To update the merge information we need to adjust the merge indicies of G
% by the number of patches in F for each level. skip(k) contains this
% number.
skip = cellfun(@(t) max(t(:)), F.mergeIdx);
mergeIdx = cell(1, treeSizeH); % Initialise.
for k = 1:(treeSizeH-1)
    sF = size(F.mergeIdx{k},2);
    sG = size(G.mergeIdx{k},2);
    if ( sF > sG )
        G.mergeIdx{k}(:,end+1) = NaN;
    elseif ( sF < sG )
        F.mergeIdx{k}(:,end+1) = NaN;
    end
    mergeIdx{k} = [F.mergeIdx{k}; G.mergeIdx{k} + skip(k)];
end
mergeIdx{end} = [1, 2]; % Final merge leaves a single patch.

if ( isnumeric(F.domain) && ~isnumeric(G.domain) )
%     F.domain = ultraSEM.rectangle(F.domain, 'makeObj');
    F.domain = util.rect2quad(F.domain);
elseif  ( ~isnumeric(F.domain) && isnumeric(G.domain) )
%     G.domain = ultraSEM.rectangle(G.domain, 'makeObj');
    G.domain = util.rect2quad(G.domain);
end

% Construct the new ULTRASEM.DOMAIN:
H = ultraSEM.Domain([F.domain ; G.domain], mergeIdx);

end
