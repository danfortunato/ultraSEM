function pass = test_updateRHS()

% Test that we can update the RHS efficiently.

tol = 1e-7;

S = ultraSEM.rectangle([-1 1 -1 1]);
T = refine(S, 2);

op = {1, 0, @(x,y) y};
rhs = -1;
bc = 0;
p = 21;

tic
S = ultraSEM(T, op, rhs, p);
S.build;
sol = S\bc;
t1 = toc;
% Obtained from CHEBOP2
pass(1) = abs(0.270753420352144 - feval(sol,0.2,0.3)) < tol;

tic
S.rhs = @(x,y) sin(pi*x.*y);
sol = S\bc;
t2 = toc;
% Obtained from CHEBOP2
pass(2) = abs(-0.017310261817778 - feval(sol,0.2,0.3)) < tol; 

pass(3) = t2 < t1/2;

%%

T = ultraSEM.triangle();

op = {1, 0.1, 1};
rhs = @(x,y) sin(pi*x.*y);
bc = 0;
p = 21;

S = ultraSEM(T, op, 0, p);
S.rhs = rhs;
sol1 = S\bc;

S2 = ultraSEM(T, op, rhs, p);
sol2 = S2\bc;

pass(4) =abs(sol1(.5,.5) - sol2(.5,.5)) < tol;

%% Test passing an ultraSEM.Sol as the RHS:

D = ultraSEM.rectangle();
D = refine(D, 1);

L = ultraSEM(D, {1,0,0}, -1, 5);
u = L\0;

L.rhs = u;
v = L\0;

L.rhs = @(x,y) u(x,y);
z = L\0;
pass(5) = normest(v-z) < tol;

end
