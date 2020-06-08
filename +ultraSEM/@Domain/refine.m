function T = refine(T, m, varargin)
%REFINE   Refine an ULTRASEM.DOMAIN.
%   REFINE(T) will divide each subdomain of T into four new equally-sized
%   pieces. The tree index information in the result is updated to reflect
%   the new subdomains, which are the first to be merged.
%
%   REFINE(T, M) will refine M times.
%
%   See also REFINEPOINT.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

if ( nargin < 2 )
    m = 1;
elseif ( size(m,2) == 2 )
    T = refinePoint(T, m, varargin{:});
    return
elseif( m == 0 ) 
    return
end

% If multiple domains, then refine each one:    
if ( numel(T) > 1 )
    for k = 1:numel(T)
        T(k) = refine(T(k), m);
    end
    return
end

% Trivial case:
if ( isempty(T.domain) )
    return
end

if ( isempty(T.mergeIdx) && length(T) > 1)
    warning('Empty tree index encountered. Building a default one.');
    T.mergeIdx = defaultIdx(T.domain);
end


if ( isnumeric(T.domain) )
    T = refineRectangle(T, m);
elseif ( isa(T.domain, 'ultraSEM.Domain') )
    for k = 1:numel(T.domain)
        T.domain(k) = refine(T.domain(k), m);
    end
else
    [T.domain, newIdx] = refine(T.domain, m);
    T.mergeIdx = [newIdx, T.mergeIdx];
end

end
