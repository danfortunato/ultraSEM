function T = rot90(T, k)
%ROT90   Rotate an ULTRASEM.DOMAIN (clockwise) by 90 degrees.
%   T = ROT90(T) rotates the ULTRASEM.DOMAIN T clockwise by 90 degrees.
%
%   T = ROT90(T, K) rotates by K*90 degrees.
%
%   See also TRANSPOSE.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

if ( nargin == 1 )
    k = 1;
end

if ( ~isinteger(k) )
    error('ULTRASEM:DOMAIN:rot90:invalid', ...
        'K must be an integer.');
end

k = mod(k, 4);
if ( k == 1 )
    T.domain = T.domain(:, [3 4 1 2]);
    T = fliplr(T);
elseif ( k == 2 )
    T = fliplr(flipud(T)); %#ok<FLUDLR>
elseif ( k == 3 )
    T.domain = T.domain(:, [3 4 1 2]);
    T = flipud(T);
else
    error('ULTRASEM:DOMAIN:rot90:invalid', ...
        'K must be an integer.');
end

end

