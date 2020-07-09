function f = times(f, g)
%.*   Pointwise multiplication for ULTRASEM.SOL.
%   F.*G multiplies F and G, where F and G may be ULTRASEM.SOL objects or
%   scalars.
%
%   See also MTIMES, COMPOSE.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

if ( ~isa(f, 'ultraSEM.Sol') )
    % Ensure F is the ULTRASEM.SOL:
    f = times(g, f);
    return
elseif ( isa(g, 'ultraSEM.Sol' ) )
    % Multiply two ULTRASEM.SOLs:
    % TODO: Check that F and G have the same domain.
    f = compose(@times, f, g);
elseif ( isnumeric(g) && isscalar(g) )
    % Multiply ULTRASEM.SOL F by scalar G:
    f.coeffs = cellfun(@(coeffs) g*coeffs, f.coeffs, 'UniformOutput', false);
else
    error('ULTRASEM:SOL:times:invalid', ...
        'F and G must be scalars or ultraSEM.Sol objects.')
end

end
