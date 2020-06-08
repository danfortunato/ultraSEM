function H = merge(varargin)
%MERGE   Merge two or more ULTRASEM.DOMAINs.
%
%   See also OLDMERGE.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

%     H = oldMerge(varargin{:}); return % Old-style merge: 

    newDomain = vertcat(varargin{:});
    if ( numel(newDomain) == 1 )
        H = newDomain;
    else
        H = ultraSEM.Domain(newDomain);
    end

    H = flatten(H);

end
