function sol = imag(sol)
%IMAG   Imaginary part of an ULTRASEM.SOL.
%
% See also REAL.

sol = compose(sol, @imag);

end
