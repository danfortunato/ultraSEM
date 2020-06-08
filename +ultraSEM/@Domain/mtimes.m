function T = mtimes(T, c)
%*   Scale an ULTRASEM.DOMAIN.
%   C*T will scale the ULTRASEM.DOMAIN by C. C must be a scalar.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

if ( ~isa(T, 'ultraSEM.Domain') )
    % Ensure T is the domain:
    T = mtimes(c, T);
    return
elseif ( isa(c, 'ultraSEM.Domain' ) )
    % We can't multiply two domains:
    error('ULTRASEM:DOMAIN:mtimes:twodomains', ...
        'Cannot multiply (* or .*) two ultraSEM.Domains.\n')
elseif ( ~isnumeric(c) )
    error('ULTRASEM:DOMAIN:mtimes:unknown', ...
        'Cannot multiply an object of type %s by an ultraSEM.Domain.', ...
        class(c));
elseif ( ~isscalar(c) )
    error('ULTRASEM:DOMAIN:mtimes:scalar', ...
        'C must be a scalar.')
end

% Scale the domain:
T.domain = c*T.domain;

end