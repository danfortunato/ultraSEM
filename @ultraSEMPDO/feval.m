function [op, isConstant] = feval(op, x, y)
%FEVAL  Evaluate the non-constant coefficients of a PDO.
%   FEVAL(OP, X, Y) evaluates the non-constant coefficients of the PDO
%   defined by the ultraSEMPDO OP.
%
% See also INITIALIZE.

% Copyright 2018 by Nick Hale and Dan Fortunato.

isConstant = true;

for k = fieldnames(op)'
    k1 = k{1};
    if ( isa(op.(k1), 'function_handle') )
        op.(k1) = diag(feval(op.(k1), x, y));
        isConstant = false;
    end
end

end
