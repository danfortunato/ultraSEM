function f = sqrt(f)
%SQRT   Square root of an ULTRASEM.SOL.
%   SQRT(F) returns the square root of the ULTRASEM.SOL F.
%
%   See also POWER, COMPOSE.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

% Empty check:
if ( isempty(f) )
    return
end

f = compose(@sqrt, f);

end
