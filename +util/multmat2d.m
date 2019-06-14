function M = multmat2d(n, A, lambda)
%MULTMAT2D   Compute the 2D N^2 x N^2 multiplication matrix for the function
%
%   f(x,y) = sum_j sum_k  A(j,k) C^{(lambda)}_j(y)C^{(lambda)}_k(x)
%
%  Alex Townsend, June 2019.

if ( nargin < 3 )
    lambda = 0;
end

if ( lambda == 0 )
    M = multmat2d_chebT(n, A);
elseif ( lambda > 0 )
    S = util.convertmat(n, 0, lambda-1);
    A = S*A*S.';
    M = multmat2d_ultraS(n, A, lambda);
else
    error('Inappropriate ultraS parameter.')
end

end

function M = multmat2d_chebT(n, A)
%MULTMAT2D_CHEBT   Compute the 2D N^2 x N^2 multiplication matrix for the function
%
%   f(x,y) = sum_j sum_k  A(j,k) T_j(y)T_k(x)
%
%  Alex Townsend, June 2019.

m = 2*n;
M = zeros(n^2, n^2);
Mold = speye(m, m);
Mx = spdiags(ones(m,2)/2,[-1 1], m, m);
Mx(2,1) = 1;
Mj = Mx;

% j = 1:
M1 = multmat1d_chebT(n, A(1,:), Mj);
M = M + kron(M1, Mold(1:n,1:n));

for j = 2:n
    M1 = multmat1d_chebT(n, A(j,:), Mx); 
    M = M + kron(M1, Mj(1:n,1:n));
    Mnew = 2*Mx*Mj - Mold;
    Mold = Mj; Mj = Mnew;
end

end 

function M = multmat2d_ultraS(n, A, lambda)
%MULTMAT2D_ULTRAS   Compute the 2D N^2 x N^2 multiplication matrix for the function
%
%   f(x,y) = sum_j sum_k  A(j,k) C^{(lambda)}_j(y)C^{(lambda)}_k(x)
%
%  Alex Townsend, June 2019.

m = 2*n;
M = zeros(n^2, n^2);
Mold = speye(m, m);
d1 = [1 2*lambda:2*lambda+m-2] ./ [1 2*((lambda+1):lambda+m-1)];
d2 = (1:m) ./ (2*(lambda:lambda + m - 1));
B = [d2' zeros(m,1) d1'];
Mx = spdiags(B, [-1 0 1], m, m);
Mj = 2*lambda*Mx;

% j = 1:
M1 = multmat1d_ultraS(n, A(1,:), Mj, lambda);
M = M + kron(M1, Mold(1:n,1:n));

for j = 2:n
    M1 = multmat1d_ultraS(n, A(j,:), Mx, lambda); 
    M = M + kron(M1, Mj(1:n,1:n));
    Mnew = (2*(j+lambda-1)/j)*Mx*Mj - ((j+2*lambda-2)/j)*Mold;
    Mold = Mj; Mj = Mnew;
end

end
