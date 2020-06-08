function pref = subsasgn(pref, ind, val)
%SUBSASGN   Subscripted assignment for ULTRASEM.PREF.
%   P.PROP = VAL, where P is an ULTRASEM.PREF object, assigns the value VAL
%   to the ULTRASEM.PREF property PROP stored in P. If PROP is not an
%   ULTRASEM.PREF property, an error will be thrown.
%
%   ULTRASEM.PREF does not support any other subscripted assignment types,
%   including '()' and '{}'.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

switch ( ind(1).type )
    case '.'
        if ( isfield(pref.prefList, ind(1).subs) )
            if ( ultraSEM.Pref.isValidPrefVal(ind(1).subs, val) )
                pref.prefList = builtin('subsasgn', pref.prefList, ind, val);
            else
                error('ULTRASEM:PREF:subsasgn:invalidPrefVal', ...
                    'Invalid preference value.');
            end
        else
            error('ULTRASEM:PREF:subsasgn:unknownPref', ...
                'Unknown preference name.');
        end
    otherwise
        error('ULTRASEM:PREF:subsasgn:badType', ...
            'Invalid subscripted assignment type.');
end

end
