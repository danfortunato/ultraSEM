function pass = test_refinePoint

pdo = {1,0,1};
rhs = -1;
n = 11;
tol = 1e-4;

% Rectangles:

D = ultraSEM.rectangle;
S1 = ultraSEM(D, pdo, rhs, n);
sol1 = S1\0;

D2 = refine(D, [1, 1]);
pass(1) = numel(D2.domain) == 3;
S2 = ultraSEM(D2, pdo, rhs, n);
sol2 = S2\0;
pass(2) = abs(sol1(.5,-.5) - sol2(.5,-.5)) < tol;

D3 = refine(D, [1, 1], 2);
pass(3) = numel(D3.domain) == 5;
S3 = ultraSEM(D3, pdo, rhs, n);
sol3 = S3\0;
pass(4) = abs(sol3(.5,-.5) - sol3(.5,-.5)) < tol;

D = refine(D, 2);
D4 = refine(D, [-1 -1 ; 1 -1 ; 1 1 ; -1 1]);
pass(5) = numel(D4.domain) == 24;
S4 = ultraSEM(D4, pdo, rhs, n);
sol4 = S4\0;
pass(6) = abs(sol1(.25,-.25) - sol2(.25,-.25)) < tol;

% Triangles:

v = [0 0 ; 1 0 ; .5 sqrt(3)/2];
T = ultraSEM.triangle(v);
S1 = ultraSEM(T, pdo, rhs, n);
sol1 = S1\0;

T2 = refine(T, v(end,:));
pass(7) = numel(T2.domain) == 5;
S2 = ultraSEM(T2, pdo, rhs, n);
sol2 = S2\0;
pass(8) = abs(sol1(.25,.25) - sol2(.25,.25)) < tol;

T3 = refine(T, v, 3);
pass(9) = numel(T3.domain) == 21;
S3 = ultraSEM(T3, pdo, rhs, 5);
sol3 = S3\0;
pass(10) = abs(sol1(.25,.25) - sol3(.25,.25)) < tol;

end

