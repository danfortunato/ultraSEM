function f = uminus(f)
%-   Unary minus for an ULTRASEM.SOL.
%   -F negates the ULTRASEM.SOL F.
%
%   G = uminus(F) is called for the syntax '-F'.
%
% See also UPLUS.

f.coeffs = cellfun(@uminus, f.coeffs, 'UniformOutput', false);

end
