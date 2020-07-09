function f = power(f, g)
%.^   Pointwise power of an ULTRASEM.SOL.
%   F.^G returns an ULTRASEM.SOL F to the scalar power G, a scalar F to the
%   ULTRASEM.SOL power G, or an ULTRASEM.SOL F to the ULTRASEM.SOL power G.
%
%   See also COMPOSE.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

if ( isempty(f) )
    return
elseif ( isempty(g) )
    f = g;
    return
else
    f = compose(@power, f, g);
end

end
