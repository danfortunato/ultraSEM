% Penrose example

T = pentaflake(2, 0, -10);
op = {{1,0,1}, {0,0}, @(x,y) 100*(1-y)};
S = ultraSEM(T, op, -1, 10);
sol = S\0;
plot(sol)
