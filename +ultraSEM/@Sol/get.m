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
% See also SUBSREF.

% Get the properties.
switch ( prop )
    case 'domain'
        val = sol.domain;
    case 'u'
        val = sol.u;
    case 'nplotpts'
        val = sol.nplotpts;
    otherwise
        error('ULTRASEM:SOL:get:propName', ...
            [prop,' is not a valid ULTRASEM.SOL property.'])
end

end
