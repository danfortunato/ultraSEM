function M = multmat1d_chebT_clenshaw( n, a, Mx )
%MULTMAT1D_CHEBT_CLENSHAW   Compute the 1D N x N multiplication matrix for the function
%
%   f(x) = sum_j a(j) T_j(x)

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
        Mnew = a(k)*I + twoMx*Ms - Mold;
        Mold = Ms; Ms = Mnew;
    end
    M = a(1)*I + Mx*Ms - Mold;
    M = M(1:n,1:n);
end

end
