function pass = test_triangle()

% Test Helmholtz on an equilateral triangle using kites.
% The RHS is chosen so that the solution is known eactly.

b = 100;
n = 10;
tol = 1e-7;

if ( ~exist('chebfun2.m', 'file') )
    pass = true;
    return
end

% Use Chebfun2 to form a suitable RHS:
u = @(x,y) exp(x.*y).*(y-sqrt(3)*x).*y.*(y+sqrt(3)*(x-1));
U = chebfun2(u, [0 1 0 1]);
F = lap(U) + b*U;
rhs = @(x,y) F(x,y);

% Set up the PDE:
op = {{1,0,1}, {0,0}, b};
T = ultraSEM.triangle([0 0; 1 0; 1/2 sqrt(3)/2]);
S = ultraSEM(T, op, rhs, n);
sol = S\0;
figure, plot(sol)

% Check the error:
err = [0 0 0];
[x,y] = getGrid(sol);
for k = 1:numel(sol.u)
    e = abs(u(x{k},y{k}) - feval(sol,x{k},y{k}));
    err(k) = norm(e, inf);
end

%%

pass = err < tol;

end
