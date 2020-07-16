function f = abs(f)
%ABS   Absolute value of an ULTRASEM.SOL.
%   ABS(F) returns the absolute value of the ULTRASEM.SOL F.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

% Empty check:
if ( isempty(f) )
    return
end

f = compose(@abs, f);

end
