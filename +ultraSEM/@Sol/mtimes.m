function f = mtimes(f, c)
%*   Scale an ULTRASEM.SOL.
%   c*F or F*c multiplies an ULTRASEM.SOL F by a scalar c.
%
%   See also TIMES.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

if ( ~isa(f, 'ultraSEM.Sol') )
    % Ensure F is the ULTRASEM.SOL:
    f = mtimes(c, f);
    return
elseif ( isa(c, 'ultraSEM.Sol' ) )
    % MTIMES should not be used to multiply two ULTRASEM.SOLs:
    error('ULTRASEM:SOL:mtimes:twosols', ...
        ['Cannot multiply two ultraSEM.Sols with ''*''. ', ...
         'Did you mean ''.*''?\n'])
elseif ( isnumeric(c) && isscalar(c) )
    % Multiply ULTRASEM.SOL F by scalar c:
    f.coeffs = cellfun(@(coeffs) c*coeffs, f.coeffs, 'UniformOutput', false);
else
    error('ULTRASEM:SOL:mtimes:invalid', 'c must be a scalar.')
end

end
