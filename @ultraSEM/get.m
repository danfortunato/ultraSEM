function val = get(S, prop)
%GET   Get properties of an ULTRASEM.
%   VAL = GET(S, PROP) returns the value of the property specified in the
%   string PROP from the ULTRASEM object S. Valid entries for the string
%   PROP are:
%
%      'DOMAIN'
%      'PATCHES'
%      'OP'
%
%   See also SUBSREF.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

% Get the properties.
switch ( prop )
    case 'domain'
        val = S.domain;
    case 'patches'
        val = S.patches;
    case 'op'
        val = S.op;
    otherwise
        error('ULTRASEM:ULTRASEM:get:propName', ...
            [prop ' is not a valid ULTRASEM property.'])
end

end
