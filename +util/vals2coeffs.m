function coeffs = vals2coeffs(values)
%VALS2COEFFS   Convert values at Chebyshev points to Chebyshev coefficients.
%
% See also COEFFS2VALS.

% Store the vals2coeffs matrices for sizes < cutoff
persistent F
cutoff = 100;

% Get the length of the input:
n = size(values, 1);

if ( n <= 1 )
    % Trivial case (constant):
    coeffs = values;
elseif ( n < cutoff )
    % Use matrix multiplication for small problems
    if ( isempty(F) )
        F = cell(cutoff, 1);
    end
    if ( isempty(F{n}) )
        x = util.chebpts(n);
        F{n} = inv( cos(acos(x) * (0:n-1)) ); % Is there a better way to do this than inv?
    end
    coeffs = F{n} * values;
else
    % Use fast transform
    coeffs = chebtech2.vals2coeffs(values); % TODO: Remove chebfun
end

end
