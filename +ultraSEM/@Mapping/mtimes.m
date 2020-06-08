function T = mtimes(T, c)
%*   Scale a mapping.
%   T*C scales the ULTRASEM.MAPPING T by the scalar C.
%
%   See also PLUS, MINUS.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

if ( isnumeric(T) ), [T, c] = deal(c, T); end
if ( ~isnumeric(c) )
    error('ULTRASEM:MAPPING:mtimes:unknown', ...
        'Cannot multiply an ultraSEM.mapping by an object of type %s.', ...
        class(c));
elseif ( ~isscalar(c) )
    error('ULTRASEM:MAPPING:mtimes:scalar', 'C must be a scalar.')
end

% Shift the domain:
for k = 1:numel(T)
    T(k).v = c*T(k).v;
end

end
