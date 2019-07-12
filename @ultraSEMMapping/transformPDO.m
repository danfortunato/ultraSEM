function [new_PDO, rhs] = transformPDO( T, PDO, rhs )
% Perform change of variables to convert PDO on mapped to domain to equivalent
% PDO on [-1,1]^2.

% Change of variables:
%
% d/dx = ds/dx d/ds + dt/dx d/dt
% d/dy = ds/dy d/ds + dt/dy d/dt
% d/dxx = (ds/dxx) d/ds + (dt/dxx) d/dt + (ds/dx)^2 d/dss +
%            (dt/dx)^2 d/dtt + 2(ds/dx)(dt/dx) d/dst
% d/dyy = (ds/dyy) d/ds + (dt/dyy) d/dt + (ds/dy)^2 d/dss +
%            (dt/dy)^2 d/dtt + 2(ds/dy)(dt/dy) d/dst
% d/dxy = (ds/dxy) d/ds + (dt/dxy) d/dt + (ds/dx)(ds/dy)d/dss +
%            (dt/dx)(dt/dy) d/dtt + ((ds/dy)(dt/dx)+(ds/dx)(dt/dy))d/dst

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEVELOPER NOTE: Manipulating function handles is expensive, so we use strings.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dsdx  = 'T.dinvT11( T.T1(s,t), T.T2(s,t) )';
dsdy  = 'T.dinvT12( T.T1(s,t), T.T2(s,t) )';
dtdx  = 'T.dinvT21( T.T1(s,t), T.T2(s,t) )';
dtdy  = 'T.dinvT22( T.T1(s,t), T.T2(s,t) )';
dsdxx = 'T.d2invT11( T.T1(s,t), T.T2(s,t) )';
dsdxy = 'T.d2invT12( T.T1(s,t), T.T2(s,t) )';
dsdyy = 'T.d2invT13( T.T1(s,t), T.T2(s,t) )';
dtdxx = 'T.d2invT21( T.T1(s,t), T.T2(s,t) )';
dtdxy = 'T.d2invT22( T.T1(s,t), T.T2(s,t) )';
dtdyy = 'T.d2invT23( T.T1(s,t), T.T2(s,t) )';

dxx = PDO.dxx;
dxy = PDO.dxy;
dyy = PDO.dyy;
dx  = PDO.dx;
dy  = PDO.dy;
b   = PDO.b;

% d/dxx term
if ( isa( dxx, 'function_handle') )
    dxx = 'dxx(T.T1(s,t), T.T2(s,t))';
else
    dxx = num2str(dxx);
end
new_dxx = ['@(s,t) ' dxx '.* ' dsdx '.^2'];
new_dxy = ['@(s,t) 2*' dxx '.* ' dsdx '.*' dtdx ''];
new_dyy = ['@(s,t) ' dxx '.*' dtdx '.^2'];
new_dx  = ['@(s,t) ' dxx '.*' dsdxx ''];
new_dy  = ['@(s,t) ' dxx '.*' dtdxx ''];

% d/dxy term
if ( isa( dxy, 'function_handle') )
    dxy = 'dxy(T.T1(s,t), T.T2(s,t))';
else
    dxy = num2str(dxy);
end
new_dxx = [new_dxx ' + ' dxy '.*' dsdx '.*' dsdy ''];
new_dxy = [new_dxy ' + 2*' dxy '.*(' dsdy '.*' dtdx '+' dsdx '.*' dtdy ')'];
new_dyy = [new_dyy ' + ' dxy '.*' dtdx '.*' dtdy ''];
new_dx  = [new_dx  ' + ' dxy '.*' dsdxy ''];
new_dy  = [new_dy  ' + ' dxy ' .*' dtdxy ''];

% d/dyy term
if ( isa( dyy, 'function_handle') )
    dyy = 'dyy(T.T1(s,t), T.T2(s,t))';
else
    dyy = num2str(dyy);
end
new_dxx = [new_dxx ' + ' dyy '.*' dsdy '.^2'];
new_dxy = [new_dxy ' + 2*' dyy '.*' dsdy '.*' dtdy ''];
new_dyy = [new_dyy ' + ' dyy '.*' dtdy '.^2'];
new_dx  = [new_dx  ' + ' dyy '.*' dsdyy ''];
new_dy  = [new_dy  ' + ' dyy '.*' dtdyy ''];

% d/dx term
if ( isa( dx, 'function_handle') )
    dx = 'dx(T.T1(s,t), T.T2(s,t))';
else
    dx = num2str(dx);
end
new_dx = [new_dx ' + ' dx '.*' dsdx ''];
new_dy = [new_dy ' + ' dx '.*' dtdx ''];

% d/dy term
if ( isa( dy, 'function_handle') )
    dy = 'dy(T.T1(s,t), T.T2(s,t))';
else
    dy = num2str(dy);
end
new_dx = [new_dx ' + ' dy '.*' dsdy ''];
new_dy = [new_dy ' + ' dy '.*' dtdy ''];

new_dxx = eval(new_dxx);
new_dxy = eval(new_dxy);
new_dyy = eval(new_dyy);
new_dx = eval(new_dx);
new_dy = eval(new_dy);

sing = @(s,t) (s/2+1/2).^2;
syms s t
e = simplify( sing(s,t) .* new_dxx(s,t) ); new_dxx = @(s,t) eval(e);
e = simplify( sing(s,t) .* new_dxy(s,t) ); new_dxy = @(s,t) eval(e);
e = simplify( sing(s,t) .* new_dyy(s,t) ); new_dyy = @(s,t) eval(e);
e = simplify( sing(s,t) .* new_dx(s,t) );  new_dx  = @(s,t) eval(e);
e = simplify( sing(s,t) .* new_dy(s,t) );  new_dy  = @(s,t) eval(e);

% zeroth term;
if ( isa(b, 'function_handle') )
    new_b = @(s,t) sing(s,t) .* b(T.T1(s,t), T.T2(s,t));
else
    new_b = @(s,t) sing(s,t) .* b;
end

new_PDO = ultraSEMPDO({new_dxx, new_dxy, new_dyy}, {new_dx, new_dy}, new_b);

% righthand side:
if ( isa(rhs, 'function_handle') )
    rhs = @(s,t) sing(s,t) .* rhs(T.T1(s,t), T.T2(s,t));
else
    rhs = @(s,t) sing(s,t) .* rhs;
end

end

function [new_PDO, rhs] = transformPDO_old( T, PDO, rhs )
% Perform change of variables to convert PDO on mapped to domain to equivalent
% PDO on [-1,1]^2.

% Change of variables:
%
% d/dx = ds/dx d/ds + dt/dx d/dt
% d/dy = ds/dy d/ds + dt/dy d/dt
% d/dxx = (ds/dxx) d/ds + (dt/dxx) d/dt + (ds/dx)^2 d/dss +
%            (dt/dx)^2 d/dtt + 2(ds/dx)(dt/dx) d/dst
% d/dyy = (ds/dyy) d/ds + (dt/dyy) d/dt + (ds/dy)^2 d/dss +
%            (dt/dy)^2 d/dtt + 2(ds/dy)(dt/dy) d/dst
% d/dxy = (ds/dxy) d/ds + (dt/dxy) d/dt + (ds/dx)(ds/dy)d/dss +
%            (dt/dx)(dt/dy) d/dtt + ((ds/dy)(dt/dx)+(ds/dx)(dt/dy))d/dst

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEVELOPER NOTE: This is the old code which uses function handles.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dsdx  = @(s,t) T.dinvT1{1}( T.T1(s,t), T.T2(s,t) );
dsdy  = @(s,t) T.dinvT1{2}( T.T1(s,t), T.T2(s,t) );
dtdx  = @(s,t) T.dinvT2{1}( T.T1(s,t), T.T2(s,t) );
dtdy  = @(s,t) T.dinvT2{2}( T.T1(s,t), T.T2(s,t) );
dsdxx = @(s,t) T.d2invT1{1}( T.T1(s,t), T.T2(s,t) );
dsdxy = @(s,t) T.d2invT1{2}( T.T1(s,t), T.T2(s,t) );
dsdyy = @(s,t) T.d2invT1{3}( T.T1(s,t), T.T2(s,t) );
dtdxx = @(s,t) T.d2invT2{1}( T.T1(s,t), T.T2(s,t) );
dtdxy = @(s,t) T.d2invT2{2}( T.T1(s,t), T.T2(s,t) );
dtdyy = @(s,t) T.d2invT2{3}( T.T1(s,t), T.T2(s,t) );

dxx = PDO.dxx;
dxy = PDO.dxy;
dyy = PDO.dyy;
dx  = PDO.dx;
dy  = PDO.dy;
b   = PDO.b;

% d/dxx term
if ( isa( dxx, 'function_handle') )
    c = @(s,t) dxx(T.T1(s,t), T.T2(s,t));
    new_dxx = @(s,t) c(s,t).*dsdx(s,t).^2;
    new_dxy = @(s,t) 2*c(s,t).*dsdx(s,t).*dtdx(s,t);
    new_dyy = @(s,t) c(s,t).*dtdx(s,t).^2;
    new_dx = @(s,t) c(s,t).*dsdxx(s,t);
    new_dy = @(s,t) c(s,t).*dtdxx(s,t);
else
    new_dxx = @(s,t) dxx.*dsdx(s,t).^2;
    new_dxy = @(s,t) 2*dxx.*dsdx(s,t).*dtdx(s,t);
    new_dyy = @(s,t) dxx.*dtdx(s,t).^2;
    new_dx = @(s,t) dxx.*dsdxx(s,t);
    new_dy = @(s,t) dxx.*dtdxx(s,t);
end

% d/dxy term
if ( isa( PDO.dxy, 'function_handle') )
    c = @(s,t) dxy(T.T1(s,t), T.T2(s,t));
    new_dxx = @(s,t) new_dxx(s,t) + c(s,t).*dsdx(s,t).*dsdy(s,t);
    new_dxy = @(s,t) new_dxy(s,t) + 2*c(s,t).*(dsdy(s,t).*dtdx(s,t)+dsdx(s,t).*dtdy(s,t));
    new_dyy = @(s,t) new_dyy(s,t) + c(s,t).*dtdx(s,t).*dtdy(s,t);
    new_dx = @(s,t) new_dx(s,t) + c(s,t).*dsdxy(s,t);
    new_dy = @(s,t) new_dy(s,t) + c(s,t).*dtdxy(s,t);
else
    new_dxx = @(s,t) new_dxx(s,t) + dxy.*dsdx(s,t).*dsdy(s,t);
    new_dxy = @(s,t) new_dxy(s,t) + 2*dxy.*(dsdy(s,t).*dtdx(s,t)+dsdx(s,t).*dtdy(s,t));
    new_dyy = @(s,t) new_dyy(s,t) + dxy.*dtdx(s,t).*dtdy(s,t);
    new_dx = @(s,t) new_dx(s,t) + dxy.*dsdxy(s,t);
    new_dy = @(s,t) new_dy(s,t) + dxy.*dtdxy(s,t);
end

% d/dyy term
if ( isa( PDO.dyy, 'function_handle') )
    dyy = @(s,t) dyy(T.T1(s,t), T.T2(s,t));
    new_dxx = @(s,t) new_dxx(s,t) + c(s,t).*dsdy(s,t).^2;
    new_dxy = @(s,t) new_dxy(s,t) + 2*c(s,t).*dsdy(s,t).*dtdy(s,t);
    new_dyy = @(s,t) new_dyy(s,t) + c(s,t).*dtdy(s,t).^2;
    new_dx = @(s,t) new_dx(s,t) + c(s,t).*dsdyy(s,t);
    new_dy = @(s,t) new_dy(s,t) + c(s,t).*dtdyy(s,t);
else
    new_dxx = @(s,t) new_dxx(s,t) + dyy.*dsdy(s,t).^2;
    new_dxy = @(s,t) new_dxy(s,t) + 2*dyy.*dsdy(s,t).*dtdy(s,t);
    new_dyy = @(s,t) new_dyy(s,t) + dyy.*dtdy(s,t).^2;
    new_dx = @(s,t) new_dx(s,t) + dyy.*dsdyy(s,t);
    new_dy = @(s,t) new_dy(s,t) + dyy.*dtdyy(s,t);
end

% d/dx term
if ( isa( dx, 'function_handle') )
    c = @(s,t) dx(T.T1(s,t), T.T2(s,t));
    new_dx = @(s,t) new_dx(s,t) + c(s,t).*dsdx(s,t);
    new_dy = @(s,t) new_dy(s,t) + c(s,t).*dtdx(s,t);
else
    new_dx = @(s,t) new_dx(s,t) + dx.*dsdx(s,t);
    new_dy = @(s,t) new_dy(s,t) + dx.*dtdx(s,t);
end

% d/dy term
if ( isa( dy, 'function_handle') )
    c = @(s,t) dy(T.T1(s,t), T.T2(s,t));
    new_dx = @(s,t) new_dx(s,t) + c(s,t).*dsdy(s,t);
    new_dy = @(s,t) new_dy(s,t) + c(s,t).*dtdy(s,t);
else
    new_dx = @(s,t) new_dx(s,t) + dy.*dsdy(s,t);
    new_dy = @(s,t) new_dy(s,t) + dy.*dtdy(s,t);
end

% zeroth term;
if ( isa(b, 'function_handle') )
    new_b = @(s,t) b(T.T1(s,t), T.T2(s,t));
else
    new_b = b;
end

new_PDO = ultraSEMPDO({new_dxx, new_dxy, new_dyy}, {new_dx, new_dy}, new_b);

% righthand side:
if ( isa(rhs, 'function_handle') )
    rhs = @(s,t) rhs(T.T1(s,t), T.T2(s,t));
end

end