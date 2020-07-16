function f = exp(f)
%EXP   Exponential of an ULTRASEM.SOL.
%   EXP(F) returns the exponential of the ULTRASEM.SOL F.
%
%   See also LOG.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

% Empty check:
if ( isempty(f) )
    return
end

f = compose(@exp, f);

end
