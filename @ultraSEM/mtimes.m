function u = mtimes(S, u)
%MTIMES   Forward application of an ULTRASEM object.
%   S*U applies the differential operator represented by the ULTRASEM S to
%   the ULTRASEM.SOL U. This is equivalent to S(U).
%
%   See also SUBSREF.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

if ( isa(S, 'ultraSEM') && isa(u, 'ultraSEM.Sol') )
    op = pdo2func(S.op);
    if ( nargin(op) == 1 )
        u = op(u);
    elseif ( nargin(op) == 3 )
        x = clone(@(x,y) x, u);
        y = clone(@(x,y) y, u);
        u = op(x, y, u);
    else
        error('ULTRASEM:ULTRASEM:mtimes:oparguments', ...
            'Unknown number of arguments to differential operator.');
    end
else
    error('ULTRASEM:ULTRASEM:mtimes:badInput', ...
        'Can only forward apply an ULTRASEM to an ULTRASEM.SOL.');
end

end
