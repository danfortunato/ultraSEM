function f = sin( f )
%SIN   Sine of an ULTRASEM.SOL.
%   SIN(F) returns the sine of the ULTRASEM.SOL F.
%
% See also ASIN, SIND.

% Empty check:
if ( isempty(f) )
    return
end

f = compose(@sin, f);

end
