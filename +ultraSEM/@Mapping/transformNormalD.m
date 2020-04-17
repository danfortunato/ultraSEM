function normal_d = transformNormalD(T, n)
%TRANSFORMNORMALD   Normal derivative operator for mapped domains.

% Chebyshev differentiation matrices in (s,t), i.e. mapped [-1,1]^2 space
I = speye( n );
D = util.diffmat( n, 1 );
S = util.convertmat( n, 0, 0 );
Dcheb = S\D;
Ds = kron( Dcheb, I );
Dt = kron( I, Dcheb );

leftIdx  = sub2ind([n n], (1:n).', ones(n,1));
rightIdx = sub2ind([n n], (1:n).', n*ones(n,1));
downIdx  = sub2ind([n n], ones(n,1), (1:n).');
upIdx    = sub2ind([n n], n*ones(n,1), (1:n).');
allIdx = [leftIdx ; rightIdx ; downIdx ; upIdx];
[X, Y] = util.chebpts2(n); % Chebyshev points and grid.
X = X(allIdx); Y = Y(allIdx);

% Evaluation matrices on the four sides
lbc = kron( (-1).^(0:n-1), I );
rbc = kron( ones(1,n), I );
dbc = kron( I, (-1).^(0:n-1) );
ubc = kron( I, ones(1,n) );

% (s,t)-differentiation matrices on the four sides
Ds_l = lbc * Ds; Dt_l = lbc * Dt;
Ds_r = rbc * Ds; Dt_r = rbc * Dt;
Ds_d = dbc * Ds; Dt_d = dbc * Dt;
Ds_u = ubc * Ds; Dt_u = ubc * Dt;

% ^^ Everything above can be precomputed for a given n ^^

% Boundary points on the four sides
[x, y] = transformGrid(T, X, Y);
xxl = x(1:n,1);           yyl = y(1:n,1);
xxr = x((n+1):2*n,1);     yyr = y((n+1):2*n,1);
xxd = x((2*n+1):3*n,1);   yyd = y((2*n+1):3*n,1);
xxu = x((3*n+1):4*n,1);   yyu = y((3*n+1):4*n,1);

% Derivatives of (s,t) with respect to (x,y)
dsdx = @(x,y) T.dinvT11(x,y);
dsdy = @(x,y) T.dinvT12(x,y);
dtdx = @(x,y) T.dinvT21(x,y);
dtdy = @(x,y) T.dinvT22(x,y);

% Get the 1D Chebyshev coefficients of d(s,t)/d(x,y) on the
% four sides through evaluation at the boundary nodes, and
% construct multiplication matrices
coeffMult = @(vals) util.multmat(n, util.vals2coeffs(vals), 0);
dsdx_l = coeffMult( dsdx(xxl,yyl) ); dsdy_l = coeffMult( dsdy(xxl,yyl) );
dtdx_l = coeffMult( dtdx(xxl,yyl) ); dtdy_l = coeffMult( dtdy(xxl,yyl) );
dsdx_r = coeffMult( dsdx(xxr,yyr) ); dsdy_r = coeffMult( dsdy(xxr,yyr) );
dtdx_r = coeffMult( dtdx(xxr,yyr) ); dtdy_r = coeffMult( dtdy(xxr,yyr) );
dsdx_d = coeffMult( dsdx(xxd,yyd) ); dsdy_d = coeffMult( dsdy(xxd,yyd) );
dtdx_d = coeffMult( dtdx(xxd,yyd) ); dtdy_d = coeffMult( dtdy(xxd,yyd) );
dsdx_u = coeffMult( dsdx(xxu,yyu) ); dsdy_u = coeffMult( dsdy(xxu,yyu) );
dtdx_u = coeffMult( dtdx(xxu,yyu) ); dtdy_u = coeffMult( dtdy(xxu,yyu) );

% (x,y)-differentiation matrices on the four sides, computed
% via the chain rule
Dx_l = dsdx_l * Ds_l + dtdx_l * Dt_l;
Dy_l = dsdy_l * Ds_l + dtdy_l * Dt_l;
Dx_r = dsdx_r * Ds_r + dtdx_r * Dt_r;
Dy_r = dsdy_r * Ds_r + dtdy_r * Dt_r;
Dx_d = dsdx_d * Ds_d + dtdx_d * Dt_d;
Dy_d = dsdy_d * Ds_d + dtdy_d * Dt_d;
Dx_u = dsdx_u * Ds_u + dtdx_u * Dt_u;
Dy_u = dsdy_u * Ds_u + dtdy_u * Dt_u;

% Convert to normal derivatives on each side
vn = normals(T);
normal_d = [ vn(1,1)*Dx_l + vn(2,1)*Dy_l ;  % "Left"
             vn(1,2)*Dx_r + vn(2,2)*Dy_r ;  % "Right"
             vn(1,3)*Dx_d + vn(2,3)*Dy_d ;  % "Down"
             vn(1,4)*Dx_u + vn(2,4)*Dy_u ]; % "Up"

end
