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

end