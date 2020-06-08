function out = subsref(pref, ind)
%SUBSREF   Subscripted reference for ULTRASEM.PREF.
%   P.PROP, where P is an ULTRASEM.PREF object, returns the value of the
%   ULTRASEM.PREF property PROP stored in P.  If PROP is not an
%   ULTRASEM.PREF property, an error will be thrown.
%
%   ULTRASEM.PREF does not support any other subscripted reference types,
%   including '()' and '{}'.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

switch ( ind(1).type )
    case '.'
        if ( isfield(pref.prefList, ind(1).subs) )
            out = pref.prefList.(ind(1).subs);
        else
            error('ULTRASEM:PREF:subsref:unknownPref', ...
                'Unknown preference name.');
        end
        if ( numel(ind) > 1 )
            out = subsref(out, ind(2:end));
        end
    otherwise
        error('ULTRASEM:PREF:subsref:badType', ...
            'Invalid subscripted reference type.');
end

end
