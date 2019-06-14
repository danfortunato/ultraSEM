function M = multmat1d_ultraS(n, a, Mx, lambda)
%MULTMAT1D_ULTRAS   Compute the 1D N x N multiplication matrix for the function
%
%   f(x,y) = sum_j a(j) C^{(lambda)}_j(x)
%
%  Alex Townsend, June 2019.

m = 2*n;
Mold = speye( m );
Ms = 2*lambda*Mx;
M = a(1)*Mold + a(2)*Ms;
for k = 2:n-1
    Mnew = (2*(k+lambda-1)/k)*Mx*Ms - ((k+2*lambda-2)/k)*Mold;
    M = M + a(k+1)*Mnew;
    Mold = Ms; Ms = Mnew;
end
M = M(1:n,1:n);

end
