function pass = test_updateRHS()

% Test that we can update the RHS efficiently.

tol = 1e-7;

S = ultraSEM.rectangle([-1 1 -1 1]);
T = refine(S, 3);

op = {1, 0, @(x,y) y};
rhs = -1;
bc = 0;
n = 21;

tic
S = ultraSEM(T, op, rhs, n);
S.build;
sol = S\bc;
t1 = toc;
err(1) = abs(0.270753420352144 - feval(sol,0.2,0.3)); % Obtained from CHEBOP2

%%

tic
S.rhs = @(x,y) sin(pi*x.*y);
sol = S\bc;
t2 = toc;
err(2) = abs(-0.017310261817778 - feval(sol,0.2,0.3)); % Obtained from CHEBOP2

%%

pass = err < tol;
pass(3) = t2 < t1/2;

[t1 , t2]

end
