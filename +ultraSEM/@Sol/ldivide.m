function f = ldivide(f, g)
%.\   Pointwise left divide for ULTRASEM.SOL.
%   F./G divides G by F, where F and G may be ULTRASEM.SOL objects or
%   scalars.
%
%   See also RDIVIDE, COMPOSE.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

f = rdivide(g, f);

end
