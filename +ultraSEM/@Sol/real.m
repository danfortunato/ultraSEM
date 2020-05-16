function f = real(f)
%REAL   Real part of and ULTRASEM.SOL.
%   REAL(F) returns the real part of the ULTRASEM.SOL object F.
%
% See also IMAJ, ABS.

f.u = cellfun(@real, f.u, 'UniformOutput', false);

end