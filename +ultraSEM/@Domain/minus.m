function T = minus(T, c)
%-   Shift an ULTRASEM.DOMAIN or setminus two ULTRASEM.DOMAINs.
%   T - C, where T is a ULTRASEM.DOMAIN and C is a scalar, is equivalent to
%   T + (-C).
%
%   T - S, where S and T are both ULTRASEM.DOMAINs, will remove common
%   patches of S and T from T.
%
%   See also PLUS.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

if ( isnumeric(T) )
    T = plus(c, (-T));
elseif ( isnumeric(c) )
    T = plus(T, (-c));
% elseif ( isa(T,'ultraSEM.Domain') && isa(c,'ultraSEM.Domain') )
%     % TODO: Allow us to cut out shapes from a domain.
%     [~, iT, ~] = intersect(T, c);
%     T = removePatch(T, iT);
else
    error('ULTRASEM:DOMAIN:plus:unknown', ...
        'Cannot subtract a %s from a %s.', class(c), class(T));
end

end