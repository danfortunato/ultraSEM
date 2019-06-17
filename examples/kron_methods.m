% Solve
%
%    lap(u) + v(x,y)*u = f on [-1,1]^2
%                    u = g on boundary
%
% by imposing BCs in four ways:
%
%   (1) Sylvester constraints (e.g. what ultraSEM & chebop2 do)
%   (2) Remove 4 "corner" DOFs
%   (3) Overdetermined system with 1D multmats
%   (4) Overdetermined system with 2D multmat
%
% Method (4) is the only one that both matches the error from (1) and
% avoids computing a low-rank representation for multmats.
%
% Dan Fortunato, June 2019.

figure(1), set(gcf, 'Position',  [0, 600, 1200, 400])

n = 50;
f = chebfun2(@(x,y) x.*sin(x.*y));
bc = 0;
%v = @(x,y) 1+0*x;
v = @(x,y) 10*sin(10*x.*y); % Nontrivial variable-coefficient
vfun = chebfun2(v);

%% True solution from chebop2
L = chebop2(@(x,y,u) diff(u,2,1) + diff(u,2,2) + vfun(x,y).*u);
L.bc = bc;
sol = solvepde(L, f, 65, 65); % Max chebop2 size

%% (1) ultraSEM solution, which imposes BCs via Sylvester constraints
tic
dom = ultraSEM.rectangle([-1 1 -1 1]);
S = ultraSEM(dom, {1, 0, vfun}, f, n);
sem = S \ bc;
toc

sem = chebfun2(sem.u{1}, 'coeffs');

figure(1), subplot(141)
plot(sem-sol), view(3), title(['ultraSEM (1): ', num2str(norm(sem-sol))])
axis square
shg

%% (2) Solve an n^2-4 x n^2-4 system and set the 4 highest coefficients to 0.
tic
ii = 1:n^2; ii([n^2-n-1 n^2-n n^2-1 n^2]) = [];
I = speye( n );
corners = [n-1 n];
lbc = kron( (-1).^(0:n-1), I );
rbc = kron( ones(1,n), I );
dbc = kron( I, (-1).^(0:n-1) );
ubc = kron( I, ones(1,n) );
% Remove the corner conditions
lbc(corners,:) = []; rbc(corners,:) = [];
dbc(corners,:) = []; ubc(corners,:) = [];
bcrows = [lbc; rbc; dbc; ubc];

% D2*X*S02.' + S02*X*D2.' + S02*My*X*(S02*Mx).' + ... = F
S02 = util.convertmat(n, 0, 1); S  = S02(1:end-2,:);
D2  = util.diffmat(n, 2);       D2 =  D2(1:end-2,:);
A = kron(S, D2) + kron(D2, S);

% Add the multiplication matrices
[C, D, R] = cdr(vfun);
C = chebcoeffs(C, n);
R = chebcoeffs(R, n);
for r = 1:length(vfun)
    Mx = ultraS.multmat( n, D(r,r) * R(:,r), 0 );
    My = ultraS.multmat( n,          C(:,r), 0 );
    SMx = S02*Mx; SMx = SMx(1:end-2,:);
    SMy = S02*My; SMy = SMy(1:end-2,:);
    A = A + kron(SMx, SMy);
end
A = [ bcrows(:,ii) ; A(:,ii) ];

BC = zeros(4*(n-2),1);
F = coeffs2(f, n, n);
F = S02*F*S02.'; F = F(1:n-2,1:n-2);
rhs = [ BC ; F(:) ];

x = A \ rhs;
x([n^2-n-1 n^2-n n^2-1 n^2]) = [0 0 0 0];
toc

X = reshape(x, n, n);
u = chebfun2(X, 'coeffs');

figure(1), subplot(142)
plot(u-sol), view(3), title(['Corner method (2): ', num2str(norm(u-sol))])
axis square
shg

%% (3) Solve an n^2+4 x n^2 overdetermined system but build multmats using low-rank slices.
tic
I = speye( n );
lbc = kron( (-1).^(0:n-1), I );
rbc = kron( ones(1,n), I );
dbc = kron( I, (-1).^(0:n-1) );
ubc = kron( I, ones(1,n) );
bcrows = [lbc; rbc; dbc; ubc];

% D2*X*S02.' + S02*X*D2.' + S02*My*X*(S02*Mx).' + ... = F
S02 = util.convertmat(n, 0, 1); S  = S02(1:end-2,:);
D2  = util.diffmat(n, 2);       D2 =  D2(1:end-2,:);
A = kron(S, D2) + kron(D2, S);

% Add the multiplication matrices
[C, D, R] = cdr(vfun);
C = chebcoeffs(C, n);
R = chebcoeffs(R, n);
for r = 1:length(vfun)
    Mx = ultraS.multmat( n, D(r,r) * R(:,r), 0 );
    My = ultraS.multmat( n,          C(:,r), 0 );
    SMx = S02*Mx; SMx = SMx(1:end-2,:);
    SMy = S02*My; SMy = SMy(1:end-2,:);
    A = A + kron(SMx, SMy);
end
A = [ bcrows ; A ];

BC = zeros(4*n,1);
F = coeffs2(f, n, n);
F = S02*F*S02.'; F = F(1:n-2,1:n-2);
rhs = [ BC ; F(:) ];

x = A \ rhs;
toc

X = reshape(x, n, n);
u = chebfun2(X, 'coeffs');

figure(1), subplot(143)
plot(u-sol), view(3), title(['Overdetermined method (3): ', num2str(norm(u-sol))])
axis square
shg

%% (4) Solve an n^2+4 x n^2 overdetermined system but build multmat directly in Kronecker form.
tic
I = speye( n );
lbc = kron( (-1).^(0:n-1), I );
rbc = kron( ones(1,n), I );
dbc = kron( I, (-1).^(0:n-1) );
ubc = kron( I, ones(1,n) );
bcrows = [lbc; rbc; dbc; ubc];

% D2*X*S02.' + S02*X*D2.' + S02*My*X*(S02*Mx).' = F
S02 = util.convertmat(n, 0, 1); S  = S02(1:end-2,:);
D2  = util.diffmat(n, 2);       D2 =  D2(1:end-2,:);

% Use the identity:
%
%   kron(Sx*Mx*Dx, Sy*My*Dy) = kron(Sx,Sy) * kron(Mx,My) * kron(Dx,Dy)
%                                           |---- M ----|
V = coeffs2(vfun, n, n);
M = util.multmat2d(n, V, 0, 0);
S02M = kron(S02,S02) * M;

% Remove 4n-4 DOFs
ii = ones(n); ii(n-1:n,:) = 0; ii = logical(kron(ii,ii));
S02M(~ii) = []; % Faster than S02M = S02M(ii)
S02M = reshape(S02M, (n-2)^2, n^2);
A = kron(S, D2) + kron(D2, S) + S02M;

A = [ bcrows ; A ];

BC = zeros(4*n,1);
F = coeffs2(f, n, n);
S02 = util.convertmat(n, 0, 1);
F = S02*F*S02.'; F = F(1:n-2,1:n-2);
rhs = [ BC ; F(:) ];

x = A \ rhs;
toc

X = reshape(x, n, n);
u = chebfun2(X, 'coeffs');

figure(1), subplot(144)
plot(u-sol), view(3), title(['Overdetermined method (4): ', num2str(norm(u-sol))])
axis square
shg
