function M = multmat1d_chebT( n, a, Mx )
%MULTMAT1D_CHEBT   Compute the 1D N x N multiplication matrix for the function
%
%   f(x) = sum_j a(j) T_j(x)

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

na = numel(a);

if ( na == 1 )
    M = a(1)*speye(n);
else
    m = 2*n;
    Mold = speye(m);
    Ms = Mx;
    M = a(1)*Mold + a(2)*Ms;
    for k = 2:na-1
        Mnew = 2*Mx*Ms - Mold;
        M = M + a(k+1)*Mnew;
        Mold = Ms; Ms = Mnew;
    end
    M = M(1:n,1:n);
end

end
