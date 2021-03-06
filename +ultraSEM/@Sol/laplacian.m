function L = laplacian(f)
%LAPLACIAN   Laplacian of an ULTRASEM.SOL.
%   LAPLACIAN(F) returns an ULTRASEM.SOL representing the Laplacian of F.
%
%   See also LAP.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

L = diff(f, 2, 2) + diff(f, 2, 1);

end
