function normal_d = transformNormalD(T, n)
%TRANSFORMNORMALD   Normal derivative operator for mapped domains.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

% Chebyshev differentiation matrices in (s,t), i.e. mapped [-1,1]^2 space
I = speye( n );
D = util.diffmat( n, 1 );
S = util.convertmat( n, 0, 0 );
Dcheb = S\D;
Ds = kron( Dcheb, I );
Dt = kron( I, Dcheb );

x = util.chebpts(n);
xxl = -ones(n,1); yyl = x;
xxr =  ones(n,1); yyr = x;
xxd = x; yyd = -ones(n,1);
xxu = x; yyu =  ones(n,1);

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

% Derivatives of (s,t) with respect to (x,y)
drdx = @(x,y) T.drdx(x,y);
drdy = @(x,y) T.drdy(x,y);
dsdx = @(x,y) T.dsdx(x,y);
dsdy = @(x,y) T.dsdy(x,y);

% Get the 1D Chebyshev coefficients of d(s,t)/d(x,y) on the
% four sides through evaluation at the boundary nodes, and
% construct multiplication matrices
coeffMult = @(vals) util.multmat(n, util.vals2coeffs(vals), 0);
drdx_l = coeffMult( drdx(xxl,yyl) ); drdy_l = coeffMult( drdy(xxl,yyl) );
dsdx_l = coeffMult( dsdx(xxl,yyl) ); dsdy_l = coeffMult( dsdy(xxl,yyl) );
drdx_r = coeffMult( drdx(xxr,yyr) ); drdy_r = coeffMult( drdy(xxr,yyr) );
dsdx_r = coeffMult( dsdx(xxr,yyr) ); dsdy_r = coeffMult( dsdy(xxr,yyr) );
drdx_d = coeffMult( drdx(xxd,yyd) ); drdy_d = coeffMult( drdy(xxd,yyd) );
dsdx_d = coeffMult( dsdx(xxd,yyd) ); dsdy_d = coeffMult( dsdy(xxd,yyd) );
drdx_u = coeffMult( drdx(xxu,yyu) ); drdy_u = coeffMult( drdy(xxu,yyu) );
dsdx_u = coeffMult( dsdx(xxu,yyu) ); dsdy_u = coeffMult( dsdy(xxu,yyu) );

% (x,y)-differentiation matrices on the four sides, computed
% via the chain rule
Dx_l = drdx_l * Ds_l + dsdx_l * Dt_l;
Dy_l = drdy_l * Ds_l + dsdy_l * Dt_l;
Dx_r = drdx_r * Ds_r + dsdx_r * Dt_r;
Dy_r = drdy_r * Ds_r + dsdy_r * Dt_r;
Dx_d = drdx_d * Ds_d + dsdx_d * Dt_d;
Dy_d = drdy_d * Ds_d + dsdy_d * Dt_d;
Dx_u = drdx_u * Ds_u + dsdx_u * Dt_u;
Dy_u = drdy_u * Ds_u + dsdy_u * Dt_u;

% Convert to normal derivatives on each side
vn = normals(T);
normal_d = [ vn(1,1)*Dx_l + vn(2,1)*Dy_l ;  % "Left"
             vn(1,2)*Dx_r + vn(2,2)*Dy_r ;  % "Right"
             vn(1,3)*Dx_d + vn(2,3)*Dy_d ;  % "Down"
             vn(1,4)*Dx_u + vn(2,4)*Dy_u ]; % "Up"

end
