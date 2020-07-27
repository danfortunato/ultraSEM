function sol = imag(sol)
%IMAG   Imaginary part of an ULTRASEM.SOL.
%
%   See also REAL.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

sol = compose(@imag, sol);

end
