function pass = test_square

tol = 1e-6;

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

%% Variable Helmholtz

D = ultraSEM.rectangle([-1 1 -1 1]);
op = {1, 0, @(x,y) sin(x)+cos(pi*y)};
rhs = @(x,y) -exp(-10*(x.^2+(y-.5).^2));
n = 21;

% N = chebop2(@(x,y,u) laplacian(u) + op{3}(x,y).*u, [-1 1 -1 1]);
% N.bc = 0;
% u = N\chebfun2(rhs);
% trueSol = feval(u, 0, 0);
trueSol = 0.042722964212428; % Obtained from CHEBOP2

S = ultraSEM(D, op, rhs, n);
sol = S\0;
u0 = feval(sol, 0, 0);
err(3) = abs(u0 - trueSol);

S = ultraSEM(refine(D, 1), op, rhs, n);
sol = S\0;
u0 = feval(sol, 0, 0);
err(4) = abs(u0 - trueSol);

%% Variable Helmholtz with advection

D = ultraSEM.rectangle([-1 1 -1 1]);
op = {1, {1,0}, @(x,y) sin(x)+cos(pi*y)};
rhs = @(x,y) -exp(-10*(x.^2+(y-.5).^2));
n = 21;

% N = chebop2(@(x,y,u) ...
%     laplacian(u) + ...
%     op{2}{1}*diffx(u) + ...
%     op{2}{2}*diffy(u) + ...
%     op{3}(x,y).*u, [-1 1 -1 1]);
% N.bc = 0;
% u = N\chebfun2(rhs);
% trueSol = feval(u, 0, 0);
trueSol = 0.040581289036653; % Obtained from CHEBOP2

S = ultraSEM(D, op, rhs, n);
sol = S\0;
u0 = feval(sol, 0, 0);
err(5) = abs(u0 - trueSol);

S = ultraSEM(refine(D, 1), op, rhs, n);
sol = S\0;
u0 = feval(sol, 0, 0);
err(6) = abs(u0 - trueSol);

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
err(7) = abs(u0 - trueSol);

S = ultraSEM(refine(D, 1), op, rhs, n);
sol = S\0;
u0 = feval(sol, 0, 0);
err(8) = abs(u0 - trueSol);

%%
pass = err < tol;

end


