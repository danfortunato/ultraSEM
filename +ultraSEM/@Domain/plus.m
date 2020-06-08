function T = plus(T, c)
%+   Shift an ULTRASEM.DOMAIN.
%   T + C will shift the ULTRASEM.DOMAIN to the right by real(c) and
%   upwards by imag(C). C must be a scalar.
%
%   See also MINUS.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

if ( ~isa(T, 'ultraSEM.Domain') )
    % Ensure T is the ULTRASEM.DOMAIN:
    T = plus(c, T);
    return
elseif ( isa(c, 'ultraSEM.Domain' ) )
    % We can't add two ULTRASEM.DOMAINs:
    error('ULTRASEM:DOMAIN:plus:twodomains', ...
        ['Cannot add (+) two ultraSEM.Domains.\n', ...
           'Did you mean to merge them with &?'])
elseif ( ~isnumeric(c) )
    error('ULTRASEM:DOMAIN:plus:unknown', ...
        'Cannot add an object of type %s to a ultraSEM.Domain.', ...
        class(c));
elseif ( ~isscalar(c) )
    error('ULTRASEM:DOMAIN:plus:scalar', ...
        'C must be a scalar.')
end

% Shift the domain:
T.domain = T.domain + c;

end
