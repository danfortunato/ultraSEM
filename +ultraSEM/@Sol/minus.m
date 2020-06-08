function h = minus(f, g)
%-   Minus for ULTRASEM.SOL.
%   F - G subtracts G from F, where F and G are ULTRASEM.SOL objects. F and
%   G must have the same domains and discretization sizes. F and G may also
%   be scalars.
%
%   See also PLUS.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

h = plus( f, (-g) );

end