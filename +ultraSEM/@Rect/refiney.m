function T = refiney(T, m)
%REFINEY   Refine a domain in the y-direction.
%   REFINEY(T) will divide each subdomain of T vertically into two new
%   equally-sized pieces. The tree index information in the result is
%   updated to reflect the new subdomains, which are the first to be
%   merged.
%
%   REFINEY(T, M) will refine M times.
%
%   See also REFINE, REFINEX.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

if ( isempty(T.domain) )
    return
end
if ( isempty(T.mergeIdx) && length(T) > 1)
    warning('Empty tree index encountered. Building a default one.');
    T.mergeIdx = defaultIdx(T.domain);
end
if ( nargin < 2 )
    m = 1;
end

for l = 1:m % Refine m times.
    nDom = size(T.domain, 1);
    dom = mat2cell(T.domain, ones(nDom, 1)); % Convert to cell.
    for k = 1:nDom
        d = dom{k}; % Subdivide this square into two new pieces:
        midPt  = mean(d(3:4));
        dom{k} = [ d(1:2), d(3), midPt ;
                   d(1:2), midPt, d(4) ];
    end
    T.domain = cell2mat(dom);              % Revert to a matrix.
    hMerge = reshape(1:2*nDom, 2, nDom).'; % New horizontal merge.
    T.mergeIdx = [hMerge, T.mergeIdx];     % Append to existing.
end

end
