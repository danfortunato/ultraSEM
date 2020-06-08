function S = subsasgn(S, index, val)
%SUBSASGN   Subscripted assignment to an ULTRASEM object.
%   S.PROP = VAL allows property assignment for S.DOMAIN and S.PATCHES.
%
%   S.RHS = F is equivalent to UPDATERHS(S, F).
%
%   S() and S{} are not supported.
%
%   See also UPDATERHS.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

switch index(1).type
    case '.'

        if ( strcmpi(index(1).subs, 'rhs') )
            S = updateRHS(S, val);
        else
            S = builtin('subsasgn', S, index, val);
        end

    otherwise
        error('ULTRASEM:ULTRASEM:subsasgn:UnexpectedType', ...
            ['??? Unexpected index.type of ' index(1).type]);
end

end
