function sol = real(sol)
%REAL   Real part of an ULTRASEM.SOL.
%
%   See also IMAG.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

sol = compose(@real, sol);

end
