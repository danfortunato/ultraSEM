function T = plus(T, c)
%+   Shift a mapping.
%   T + C will shift the ULTRASEM.MAPPING T to the right by real(C) and
%   upwards by imag(C). C must be a scalar.
%
%   See also MINUS, MTIMES.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

if ( ~isnumeric(c) )
    error('ULTRASEM:MAPPING:PLUS:unknown', ...
        'Cannot add an object of type %s to a ultraSEM.Domain.', ...
        class(c));
elseif ( ~isscalar(c) )
    error('ULTRASEM:MAPPING:PLUS:scalar', ...
        'C must be a scalar.')
end

% Shift the domain:
for k = 1:numel(T)
    T(k).v(:,1) = T(k).v(:,1) + real(c);
    T(k).v(:,2) = T(k).v(:,2) + imag(c);
end

end
