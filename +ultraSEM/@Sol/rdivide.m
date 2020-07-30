function f = rdivide(f, g)
%./   Pointwise right divide for ULTRASEM.SOL.
%   F./G divides F by G, where F and G may be ULTRASEM.SOL objects or
%   scalars.
%
%   See also LDIVIDE, COMPOSE.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

% If either f or g are empty then return an empty ULTRASEM.SOL object.
if ( isempty(f) )
    return
elseif ( isempty(g) )
    f = g;
    return
end

if ( isa(f, 'ultraSEM.Sol') && isa(g, 'ultraSEM.Sol') )
    % Divide two ULTRASEM.SOLs:
    % TODO: Check that F and G have the same domain.
    [ss, wzero] = singleSignTest(g);
    if ( ss && ~wzero )
        f = compose(@rdivide, f, g);
    else
        error('ULTRASEM:SOL:rdivide:zero', ...
              'Attempting to invert an ultraSEM.Sol with a root.');
    end
elseif ( isa(f, 'ultraSEM.Sol') && isnumeric(g) && isscalar(g) )
    % Divide ULTRASEM.SOL F by scalar G:
    f.coeffs = cellfun(@(coeffs) (1./g).*coeffs, f.coeffs, 'UniformOutput', false);
elseif ( isnumeric(f) && isscalar(f) && isa(g, 'ultraSEM.Sol') )
    % Divide scalar F by ULTRASEM.SOL G:
    [ss, wzero] = singleSignTest(g);
    if ( ss && ~wzero )
        f = compose(@rdivide, f, g);
    else
        error('ULTRASEM:SOL:rdivide:zero', ...
              'Attempting to invert an ultraSEM.Sol with a root.');
    end
else
    error('ULTRASEM:SOL:rdivide:invalid', ...
        'F and G must be scalars or ultraSEM.Sol objects.')
end

end
