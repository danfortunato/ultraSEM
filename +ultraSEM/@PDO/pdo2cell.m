function C = pdo2cell(L)
%PDO2CELL   Convert a PDO to a cell array of coefficients.
%   C = PDO2CELL(L) extracts the coefficients of the PDO L and places them
%   in a cell array C = {{dxx, dxy, dyy}, {dx, dy}, b}.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

C = cell(1, 3);
C{1} = {L.dxx, L.dxy, L.dyy};
C{2} = {L.dx, L.dy};
C{3} = L.b;

end
