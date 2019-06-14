pdo = {{1, 0, 1}, {0, 0}, 100};
U = ultraSEM.alphabet('u');
S = ultraSEM.alphabet('S')+3;
E = ultraSEM.alphabet('E')+6;
M = ultraSEM.alphabet('M')+9;
D = U & S & E & M;
op = ultraSEM(D, pdo, -1, 51 );
sol = op \ 0;
plot(sol)