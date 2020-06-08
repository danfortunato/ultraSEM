function [op, isConstant] = feval(op, x, y)
%FEVAL  Evaluate the non-constant coefficients of a PDO.
%   FEVAL(OP, X, Y) evaluates the non-constant coefficients of the PDO
%   defined by the ULTRASEM.PDO OP.
%
%   See also INITIALIZE.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

isConstant = true;

for k = fieldnames(op)'
    k1 = k{1};
    if ( isa(op.(k1), 'function_handle') )
        op.(k1) = diag(feval(op.(k1), x, y));
        isConstant = false;
    end
end

end
