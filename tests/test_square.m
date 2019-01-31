function pass = test_square

tol = 1e-5;

%% Constant Helmholtz

D = ultraSEM.rectangle([-1 1 -1 1]);
op = {1, 0, 1};
rhs = -1;
n = 21;

% N = chebop2(@(u) laplacian(u) + u, [-1 1 -1 1]);
% N.bc = 0;
% u = N\rhs;
% trueSol = feval(u, 0, 0);
trueSol = 0.376530496242783; % Obtained from CHEBOP2

S = ultraSEM(D, op, rhs, n);
sol = S\0;
u0 = feval(sol, 0, 0);
err(1) = abs(u0 - trueSol);

S = ultraSEM(refine(D, 1), op, rhs, n);
sol = S\0;
u0 = feval(sol, 0, 0);
err(2) = abs(u0 - trueSol);

%% Laplacian with cross derivatives

D = ultraSEM.rectangle([-1 1 -1 1]);
op = {{1,.5,1}, {0,0}, 0};
rhs = @(x,y) 1+0*x;
n = 41;

% N = chebop2(@(u) ...
%     laplacian(u) + ...
%     op{1}{2}*diffx(diffy(u)), [-1 1 -1 1]);
% N.bc = 0;
% u = solvepde(N,chebfun2(rhs),100,100);
% trueSol = feval(u, 0, 0);
trueSol = -0.298398200773298; % Obtained from CHEBOP2

S = ultraSEM(D, op, rhs, n);
sol = S\0;
u0 = feval(sol, 0, 0);
err(3) = abs(u0 - trueSol);

S = ultraSEM(refine(D, 1), op, rhs, n);
sol = S\0;
u0 = feval(sol, 0, 0);
err(4) = abs(u0 - trueSol);

%%
pass = err < tol;

end


