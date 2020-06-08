function M = multmat1d_ultraS_clenshaw(n, a, Mx, lambda)
%MULTMAT1D_ULTRAS   Compute the 1D N x N multiplication matrix for the function
%
%   f(x,y) = sum_j a(j) C^{(lambda)}_j(x)

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

na = numel(a);

if ( na == 1 )
    M = a(1)*speye(n);
else
    m = 2*n;
    Mold = 0*Mx;
    Ms = 0*Mx;
    twoMx = 2*Mx;
    I = speye(m);
    for k = na:-1:2
        Mnew = a(k)*I + (k+lambda-1)/(k)*twoMx*Ms - (k+2*lambda-1)/(k+1)*Mold;
        Mold = Ms; Ms = Mnew;
    end
    M = a(1)*I + lambda*twoMx*Ms - lambda*Mold;
    M = M(1:n,1:n);
end

end

