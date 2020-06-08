function val = get(sol, prop)
%GET   Get properties of an ULTRASEM.SOL.
%   VAL = GET(F, PROP) returns the value of the property specified in the
%   string PROP from the ULTRASEM.SOL object SOL. Valid entries for the
%   string PROP are:
%
%      'DOMAIN'
%      'U'
%      'NPLOTPTS'
%
%   See also SUBSREF.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

% Get the properties.
switch ( prop )
    case 'domain'
        val = sol.domain;
    case 'coeffs'
        val = sol.coeffs;
    case 'nplotpts'
        val = sol.nplotpts;
    otherwise
        error('ULTRASEM:SOL:get:propName', ...
            [prop,' is not a valid ULTRASEM.SOL property.'])
end

end
