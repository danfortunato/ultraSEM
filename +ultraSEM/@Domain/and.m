function C = and(A, B)
%&   Merge two ULTRASEM.DOMAINs.
%   C = A & B will merge the two ULTRASEM.DOMAINs A and B. See MERGE() for
%   further details.
%
%   See also MERGE.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

if ( isempty(A) )
    C = B;
    return
elseif ( isempty(B) )
    C = A;
    return
end

C = merge(A, B); % This method is simply a wrapper for MERGE().

end
