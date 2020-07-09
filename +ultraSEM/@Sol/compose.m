function f = compose(op, f, g)
%COMPOSE   Compose a function handle with an ULTRASEM.SOL.
%   COMPOSE(OP, F) returns an ULTRASEM.SOL representing OP(F), where OP is
%   a function handle and F is an ULTRASEM.SOL. The function handle OP is
%   applied to F in value space.
%
%   COMPOSE(OP, F, G) returns an ULTRASEM.SOL representing OP(F, G), where
%   OP is a function handle and at least one of F and G is an ULTRASEM.SOL.
%   The function handle OP is applied to F and G in value space.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

if ( nargin == 2 )
    ff = coeffs2vals(f.coeffs);
    ff = cellfun(op, ff, 'UniformOutput', false);
elseif ( isa(f, 'ultraSEM.Sol') && isa(g, 'ultraSEM.Sol') )
    % F and G are both ULTRASEM.SOLs:
    ff = coeffs2vals(f.coeffs);
    gg = coeffs2vals(g.coeffs);
    ff = cellfun(op, ff, gg, 'UniformOutput', false);
elseif ( isa(f, 'ultraSEM.Sol') )
    % G is not an ULTRASEM.SOL:
    ff = coeffs2vals(f.coeffs);
    ff = cellfun(@(x) op(x,g), ff, 'UniformOutput', false);
elseif ( isa(g, 'ultraSEM.Sol') )
    % F is not an ULTRASEM.SOL:
    gg = coeffs2vals(g.coeffs);
    ff = cellfun(@(x) op(f,x), gg, 'UniformOutput', false);
    f = g;
end

f.coeffs = vals2coeffs(ff);

end
