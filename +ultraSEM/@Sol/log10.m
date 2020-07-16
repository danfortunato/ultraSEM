function f = log10(f)
%LOG10   Base 10 logarithm of an ULTRASEM.SOL.
%   LOG10(F) returns the base 10 logarithm of the ULTRASEM.SOL F. If F has
%   any roots in its domain, then the representation is likely to be
%   inaccurate.
%
%   See also LOG, LOG2, EXP.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

% Empty check:
if ( isempty(f) )
    return
end

f = compose(@log10, f);

end
