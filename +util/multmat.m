function M = multmat(n, a, lambda)
%MULTMAT   Compute the 1D N x N multiplication matrix for the function
%
%   f(x) = sum_j  a(j) C^{(lambda)}_j(x)
%
%  Alex Townsend, June 2019.

if ( nargin < 3 )
    lambda = 0;
end

m = 2*n;
if ( lambda == 0 )
    Mx = spdiags(ones(m,2)/2, [-1 1], m, m);
    Mx(2,1) = 1;
    M = util.multmat1d_chebT(n, a, Mx);
elseif ( lambda > 0 )
    a = util.convertmat(n, 0, lambda-1) * a;
    d1 = [1 2*lambda:2*lambda+m-2] ./ [1 2*((lambda+1):lambda+m-1)];
    d2 = (1:m) ./ (2*(lambda:lambda + m - 1));
    B = [d2' zeros(m,1) d1'];
    Mx = spdiags(B, [-1 0 1], m, m);
    M = util.multmat1d_ultraS(n, a, Mx, lambda);
else
    error('Inappropriate ultraS parameter.')
end

end
