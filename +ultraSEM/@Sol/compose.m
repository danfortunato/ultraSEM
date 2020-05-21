function f = compose(op, f, g)
%COMPOSE   Compose a function handle with an ULTRASEM.SOL.
%   COMPOSE(OP, F) returns an ULTRASEM.SOL representing OP(F), where OP is
%   a function handle and F is an ULTRASEM.SOL. The function handle OP is
%   applied to F in value space.
%
%   COMPOSE(OP, F, G) returns an ULTRASEM.SOL representing OP(F, G), where
%   OP is a function handle and F and G are ULTRASEM.SOLs. The function
%   handle OP is applied to F and G in value space.

ff = coeffs2vals(f.u);

if ( nargin == 2 )
    ff = cellfun(op, ff, 'UniformOutput', false);
else
    gg = coeffs2vals(g.u);
    ff = cellfun(op, ff, gg, 'UniformOutput', false);
end

f.u = vals2coeffs(ff);

end
