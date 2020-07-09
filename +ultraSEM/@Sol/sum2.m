function I = sum2(f)
%SUM2   Double integral of an ULTRASEM.SOL.
%   I = SUM2(F) returns the double integral of the ULTRASEM.SOL F over its
%   domain.
%
%   I = SUM2(F, 'all') returns an array of double integrals over each patch
%   of F.
%
% See also INTEGRAL2.

I = integral2(f);

end
