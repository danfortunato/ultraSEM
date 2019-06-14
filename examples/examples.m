%% Poisson on a square

n = 40;
pdo = {1, 0, 0};
rhs = -1;
bc = 0;
dom = ultraSEM.rectangle([-1 1 -1 1]);
exampleplot(dom)

op = ultraSEM(dom, pdo, rhs, n);
sol = op \ bc;
exampleplot(sol)

%% Cross derivatives

n = 100;
pdo = {{1,2,1}, 0, 0};
rhs = -1;
bc = 0;
dom = ultraSEM.rectangle([-1 0 -1 1]) & ...
      ultraSEM.rectangle([ 0 1 -1 1]);
exampleplot(dom) 

op = ultraSEM(dom, pdo, rhs, n);
sol = op \ bc;

exampleplot(sol)

%% Letters

n = 50;
pdo = {{1,0,1}, 0, 100+10i};
rhs = @(x,y) cos(x.*y)+y.^2;
bc = @(x,y) 10*sin(x+y);

S = ultraSEM.alphabet('S') + 3i;
U = ultraSEM.alphabet('U') + 3;
SU = S & U;
SU = refine(SU, 1);
exampleplot(SU)

op = ultraSEM(SU, pdo, rhs, n);
sol = op \ bc;
exampleplot(sol)

%% ...

n = 60;
v = randnfun2(1, [-10 10 4 22]) + 20;
pdo = {{1,0,1}, 0, 20*v};
rhs = -1;
bc = 0;

H = ultraSEM.alphabet('H');
E = ultraSEM.alphabet('E');
L = ultraSEM.alphabet('L');
O = ultraSEM.alphabet('O');
HELLO = (H + 6i-6) & (E + 3i-3) & (L+1i+1) & (L-2i+4) & (O-4i+8);
HELLO = HELLO + 10i;
exampleplot(HELLO)

op = ultraSEM(HELLO, pdo, rhs, n);
sol = op \ bc;
exampleplot(sol)

%% Weird merges

clf
set(gcf, 'Position',  [0, 600, 500, 700])

pdo = {1, 0, 0};
rhs = -1;
bc = 0;

%%% Merge an O-shape and a square to form a larger square

o = ultraSEM.alphabet('o');
s = ultraSEM.rectangle([0 1 0 1]);
dom = o & s;

subplot(321)
plot(o), hold on, plot(s,'b'), hold off
axis square tight

op = ultraSEM(dom, pdo, rhs);
sol = op \ bc;
subplot(322), plot(sol), axis square tight

%%% Merge an L-shape and a square to form a larger square

l = ultraSEM.alphabet('l');
s = ultraSEM.rectangle([0 1 0 1], 1);
dom = l & s;

subplot(323)
plot(l), hold on, plot(s,'b'), hold off
axis square tight

op = ultraSEM(dom, pdo, rhs);
sol = op \ bc;
subplot(324), plot(sol), axis square tight

%%% Merge pairs anti-diagonally

S14 = ultraSEM.rectangle([-1 0 0 1]) & ultraSEM.rectangle([ 0 1 -1 0]);
S23 = ultraSEM.rectangle([ 0 1 0 1]) & ultraSEM.rectangle([-1 0 -1 0]);
dom = S14 & S23;

subplot(325)
plot(S14), hold on, plot(S23,'b'), hold off
axis square tight

op = ultraSEM(dom, pdo, rhs);
sol = op \ bc;
subplot(326), plot(sol), axis square tight
shg

%% Kite

n = 30;
pdo = {0.3, 10, 0};
rhs = -1;
bc = 0;
dom = ultraSEM.kite([3 2; 1 1; 0 0; 2 0]);
exampleplot(dom)

op = ultraSEM(dom, pdo, rhs, n);
sol = op \ bc;
exampleplot(sol)

%% ...

n = 30;
pdo = {1, 0, 0};
rhs = @(x,y) 10*x+y;
bc = 0;

K1 = ultraSEM.kite([ 0 0 ;  1/2 .1 ;  .6 .8 ; 0 sqrt(3)/6]);
K2 = ultraSEM.kite([ 0 0 ; -1/2 .2 ; -.6 .5 ; 0 sqrt(3)/6]); % anticlockwise
K3 = ultraSEM.kite([ 0 sqrt(3)/6 ; .6 .8 ; .5 1 ; -.6 .5]);
K4 = ultraSEM.kite([-1/2 .2 ; -.6 0 ; -.2 -.2 ; 0 0 ]);
dom = K1 & K2 & K3 & K4;
exampleplot(dom)

op = ultraSEM(dom, pdo, rhs, n);
sol = op \ bc;
exampleplot(sol)

%% Triangle

b = 100;
pdo = {1, 0, b};
bc = 0;
dom = ultraSEM.triangle([0 0; 1 0; 1/2 sqrt(3)/2]);
exampleplot(dom)

% Exact solution
u = @(x,y) exp(x.*y).*(y-sqrt(3)*x).*y.*(y+sqrt(3)*(x-1));
U = chebfun2(u, [0 1 0 1]);
F = lap(U) + b*U;
rhs = @(x,y) F(x,y);

op = ultraSEM(dom, pdo, rhs);
sol = op \ bc;
exampleplot(sol)

err = 0;
for k = 1:numel(sol.x)
    x = sol.x{k};
    y = sol.y{k};
    e = abs(u(x,y) - feval(sol,x,y));
    err = max(err, norm(e, inf));
end
err

%% Polygon

n = 50;
pdo = {1, 0, 1000};
rhs = -1;
bc = 0;
pentagon = nsidedpoly(5);
dom = ultraSEM.polygon(pentagon.Vertices);
exampleplot(dom)

op = ultraSEM(dom, pdo, rhs, n);
sol = op \ bc;
exampleplot(sol)

%% Penrose snowflake

n = 20;
pdo = {{1,0,1}, {0,0}, @(x,y) 30*(1-y)};
rhs = -1;
bc = 0;

T = pentaflake(2, 0, -5);       % Construct Penrose snowflake shape
T = rmholes(T);                 % Remove holes
tri = triangulation(T);         % Triangulate
dom = ultraSEM.trimesh(tri);    % Convert to ultraSEMDomain
exampleplot(dom)

tic
op = ultraSEM(dom, pdo, rhs, n); % Construct solution operator on mesh
sol = op\bc;                     % Solve using BCs
toc

tic
exampleplot(sol)
toc
