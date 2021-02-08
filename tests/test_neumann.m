function pass = test_neumann()

tol = 1e-8;

n = 20;
rhs = -1;
pdo = {1, 0, 0};
dom = ultraSEM.rectangle([-2 2 -1 1], 1, 2);
S = ultraSEM(dom, pdo, rhs, n);

% The ordering of the BCs should match the ordering of S.patches{1}.edges
bc = ultraSEM.BC;
bc(1:6) = ultraSEM.BC.dirichlet(0);
bc(4)   = ultraSEM.BC.neumann(0);

u = S \ bc;
dudx = diffx(u);

pass(1) = abs(u(2,0) - 0.498072714) < tol;
pass(2) = abs(dudx(2,0)) < tol;


%%

% Compare [     ] with [ ] and [    / using Neumann BCs.

%         [     ]      [ ]     [  /
%         [     ]              |/

tol = 1e-8;

n = 20;
rhs = -1;
pdo = {1, 0, 0};
dom = ultraSEM.rectangle([-1 1 -1 1], 1, 2);
S1 = ultraSEM(dom, pdo, rhs, n);
sol1 = S1\0;

dom = ultraSEM.rectangle([0 1 0 1]);
S2 = ultraSEM(dom, pdo, rhs, n);
bc = ultraSEM.BC;
bc(1:4) = ultraSEM.BC.dirichlet(0);
bc([1 3]) = ultraSEM.BC.neumann(0);
sol2 = S2\bc;
pass(3) = abs(sol1(.5,.5) - sol2(.5,.5)) < tol;

dom = ultraSEM.triangle([-1 -1 ; -1 1 ; 1 1]);
% dom = refine(dom);
S3 = ultraSEM(dom, pdo, rhs, n);
bc = ultraSEM.BC;
bc(1:6) = ultraSEM.BC.dirichlet(0);
bc([2 4]) = ultraSEM.BC.neumann(0);
sol3 = S3\bc;
pass(4) = abs(sol1(-.5,.5) - sol3(-.5,.5)) < 100*tol;

end