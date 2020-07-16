function f = log2(f)
%LOG2   Base 2 logarithm of an ULTRASEM.SOL.
%   LOG2(F) returns the base 2 logarithm of the ULTRASEM.SOL F. If F has
%   any roots in its domain, then the representation is likely to be
%   inaccurate.
%
%   See also LOG, LOG10, EXP.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

% Empty check:
if ( isempty(f) )
    return
end

f = compose(@log2, f);

end
