function sol = real(sol)
%REAL   Real part of an ULTRASEM.SOL.
%
% See also IMAG.

sol = compose(sol, @real);

end
