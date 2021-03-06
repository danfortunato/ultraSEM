function solf = clone(f, sol)
%CLONE   Clone an ULTRASEM.SOL.
%   CLONE(F, SOL) returns an ULTRASEM.SOL with the same domain and
%   discretization size as the ULTRASEM.SOL SOL, but with coefficients that
%   represent the function F, which may be a constant, a function handle,
%   or another ULTRASEM.SOL.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

if ( ~isa(sol, 'ultraSEM.Sol') )
    error('ULTRASEM:SOL:clone:invalid', ...
        'Second argument must be an ULTRASEM.SOL.');
end

% Copy the object.
solf = sol;

% Reinitialize coefficients with the provided function.
if ( isa(f, 'function_handle') || isa(f, 'ultraSEM.Sol') )
    for k = 1:length(sol)
        [xx, yy] = getGrid(sol, k);
        vals = feval(f, xx, yy);
        solf.coeffs{k} = util.vals2coeffs( util.vals2coeffs( vals ).' ).';
    end
else
    for k = 1:length(sol)
        solf.coeffs{k} = 0*sol.coeffs{k};
        solf.coeffs{k}(1,1) = f;
    end
end

end
