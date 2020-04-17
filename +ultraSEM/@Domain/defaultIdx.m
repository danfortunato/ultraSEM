function mergeIdx = defaultIdx(dom)
%DEFAULTIDX   Make a default (naive) merge index for a domain.
%   MERGEIDX = DEFAULTIDX(DOM) automatically constructs a merge index
%   MERGEIDX for the domain DOM. The merge index will take the form
%
%       MERGEIDX = {[1 2 ; 3 4 ; ... ; N-1 N],
%                   [1 2 ; ... ; N/2-1 N/2],
%                    ...,
%                   [1 2]},
%
%   where N is the number of patches in DOM. NaNs are introduced where
%   there are an odd number of patches on a level.
%
%   MERGEIDX = DEFAULTIDX(N) is equivalent.

% Parse input:
if ( isa(dom, 'scalar') )
    np = dom;             % Number of patches.
else
    np = size(dom, 1);    % Number of patches.
end
nl = ceil(log2(np));      % Number of layers.
mergeIdx = cell(1, nl);   % Initialise.

for l = 1:nl
    if ( mod( np, 2 ) )
        ii = [1:np, NaN]; % Pad with a NaN.
    else
        ii = 1:np;
    end
    np = ceil( np/2 );    % Number of patches is halved each layer.
    mergeIdx{1,l} = reshape(ii, 2, np).';
end

end
