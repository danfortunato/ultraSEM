function M = multmat(n, a, lambda)
%MULTMAT   Compute the 1D N x N multiplication matrix for the function
%
%   f(x) = sum_j  a(j) C^{(lambda)}_j(x)
%
%  Alex Townsend, June 2019.

if ( nargin < 3 )
    lambda = 0;
end

% Chop coefficients less than TOL
tol = eps;

% Store NSTORE multmats of size NMAX x NMAX for C{(0:2)}_(0:NSTORE-1)(x).
persistent MStore nmax nstore
nmax_limit   = 100;
nstore_limit = 100;
recompute = isempty(MStore) || ...
            isempty(nmax)   || ...
            isempty(nstore) || ...
            (nmax < n && nmax <= nmax_limit) || ...
            (nstore < n && nstore <= nstore_limit);

if ( recompute )
    nmax   = min(n, nmax_limit);   % Store multmats of size nmax x NMAX
    nstore = min(n, nstore_limit); % Store the first NSTORE multmats
    MStore = cell(nstore, 3);
    m = 2*nmax;

    Mx = spdiags(ones(m,2)/2, [-1 1], m, m);
    Mx(2,1) = 1;
    for k = 1:nstore
        MStore{k,1} = util.multmat1d_chebT_clenshaw(nmax, [zeros(k-1,1); 1], Mx);
    end

    lam = 1;
    d1 = [1 2*lam:2*lam+m-2] ./ [1 2*((lam+1):lam+m-1)];
    d2 = (1:m) ./ (2*(lam:lam+m-1));
    B = [d2' zeros(m,1) d1'];
    Mx = spdiags(B, [-1 0 1], m, m);
    for k = 1:nstore
        MStore{k,2} = util.multmat1d_ultraS_clenshaw(nmax, [zeros(k-1,1); 1], Mx, lam);
    end

    lam = 2;
    d1 = [1 2*lam:2*lam+m-2] ./ [1 2*((lam+1):lam+m-1)];
    d2 = (1:m) ./ (2*(lam:lam+m-1));
    B = [d2' zeros(m,1) d1'];
    Mx = spdiags(B, [-1 0 1], m, m);
    for k = 1:nstore
        MStore{k,3} = util.multmat1d_ultraS_clenshaw(nmax, [zeros(k-1,1); 1], Mx, lam);
    end
end

if ( lambda > 0 )
    a = util.convertmat(n, 0, lambda-1) * a;
end
a = chopCoeffs(a, tol);

if ( isempty(a) )
    M = sparse(n);
elseif ( n <= nmax && numel(a) <= nstore && lambda <= 2 )
    idx = find(abs(a) > tol);
    M = a(idx(end))*MStore{idx(end),lambda+1}(1:n,1:n);
    for k = numel(idx)-1:-1:1
        M = M + a(idx(k))*MStore{idx(k),lambda+1}(1:n,1:n);
    end
else
    m = 2*n;
    if ( lambda == 0 )
        Mx = spdiags(ones(m,2)/2, [-1 1], m, m);
        Mx(2,1) = 1;
        M = util.multmat1d_chebT_clenshaw(n, a, Mx);
    elseif ( lambda > 0 )
        d1 = [1 2*lambda:2*lambda+m-2] ./ [1 2*((lambda+1):lambda+m-1)];
        d2 = (1:m) ./ (2*(lambda:lambda+m-1));
        B = [d2' zeros(m,1) d1'];
        Mx = spdiags(B, [-1 0 1], m, m);
        M = util.multmat1d_ultraS_clenshaw(n, a, Mx, lambda);
    else
        error('Inappropriate ultraS parameter.')
    end
end

end

function a = chopCoeffs(a, tol)
%CHOPCOEFFS   Chop off trailing coefficients that are less than TOL.

    if ( nargin == 1 )
        tol = eps;
    end

    if ( ~isvector(a) )
        error('Input must be a vector of coefficients.');
    end

    na = find(abs(a) > tol, 1, 'last');

    if ( isempty(na) )
        a = [];
    else
        a = a(1:na);
    end

end
