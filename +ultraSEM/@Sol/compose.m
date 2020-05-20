function sol = compose(sol, op)
%COMPOSE   Compose a function handle with an ULTRASEM.SOL.
%   COMPOSE(SOL, OP) returns an ULTRASEM.SOL representing OP(SOL), where
%   SOL is also an ULTRASEM.SOL and OP is a function handle. The function
%   handle OP is applied to SOL in value space.

vals = coeffs2vals(sol.u);
vals = cellfun(op, vals, 'UniformOutput', false);
sol.u = vals2coeffs(vals);

end
