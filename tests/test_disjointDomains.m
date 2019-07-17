function pass = test_disjointDomains

tol = 1e-10;
p = 5;

R = ultraSEM.rectangle;
T = ultraSEM.triangle;
T2 = ultraSEM.triangle([1 1 ; 2 1 ; 1.5 1-sqrt(3)/2]);

D = .5*R & (.5*R+2);
op = ultraSEM(D, {1,0,0}, -1, p);
u = op\0;
op1 = ultraSEM(.5*R, {1,0,0}, -1, p);
u1 = op1\0;
op2 = ultraSEM((.5*R+2), {1,0,0}, -1, p);
u2 = op2\0;
pass(1) = norm([u(0,0) - u1(0,0) ; u(2,.1) - u2(2,.1)], inf) < tol;

D = T & T2;
op = ultraSEM(D, {1,0,0}, -1, p);
u = op\0;
op1 = ultraSEM(T, {1,0,0}, -1, p);
u1 = op1\0;
op2 = ultraSEM(T2, {1,0,0}, -1, p);
u2 = op2\0;
pass(2) = norm([u(.2,.2) - u1(0.2,0.2) ; u(1.5,0.5) - u2(1.5,.5)], inf) < tol;

D1 = T;
D2 = (1/2.5*R-.5+.4i);
D = D1 & D2;
op = ultraSEM(D, {1,0,0}, -1, p);
u = op\0;
op1 = ultraSEM(D1, {1,0,0}, -1, p);
u1 = op1\0;
op2 = ultraSEM(D2, {1,0,0}, -1, p);
u2 = op2\0;
pass(3) = norm([u(.5,.5) - u1(.5,.5) ; u(-.5,.5) - u2(-.5,.5)], inf) < tol;

end

