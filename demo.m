close all
clc

S = ultraSEM.alphabet('S') + 2i;
U = ultraSEM.alphabet('U') + 2;
SU = S & U;
SU = refine(SU, 1);
pdo = {{1, 0, 1}, {0, 0}, 10};
rhs = 1;
bc = 0; %bc = @(x,y) x+y;
n = 30;
op = ultraSEM(SU, pdo, rhs, n);
tic
build(op)
toc
sol = op \ bc;
plot(sol), colorbar, axis off
