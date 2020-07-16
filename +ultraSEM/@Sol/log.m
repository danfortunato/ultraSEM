function f = log(f)
%LOG   Natural logarithm of an ULTRASEM.SOL.
%   LOG(F) returns the natural logarithm of the ULTRASEM.SOL F. If F has
%   any roots in its domain, then the representation is likely to be
%   inaccurate.
%
%   See also LOG2, LOG10, EXP.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

% Empty check:
if ( isempty(f) )
    return
end

f = compose(@log, f);

end
