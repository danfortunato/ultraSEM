function w = quadwts(n)
%QUADWTS   Quadrature weights for 2nd-kind Chebyshev points.
%   QUADWTS(N) returns the N weights for Clenshaw-Curtis quadrature on
%   2nd-kind Chebyshev points.
%
% See also UTIL.CHEBPTS.

if ( n == 0 )
    w = [];
elseif ( n == 1 )
    w = 2;
else
    c = 2./[1, 1-(2:2:(n-1)).^2];  % Exact integrals of T_k (even)
    c = [c, c(floor(n/2):-1:2)];   % Mirror for DCT via FFT
    w = ifft(c);                   % Interior weights
    w([1,n]) = w(1)/2;             % Boundary weights
end

end
