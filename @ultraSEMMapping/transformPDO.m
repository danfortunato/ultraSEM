function [op_rs, rhs_rs] = transformPDO( dom, op_xy, rhs_xy )
%TRANSFORMPDO Perform change of variables to convert PDO on mapped
%(x,y)-domain to equivalent PDO on reference (r,s)-domain [-1,1]^2.

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

det   = 'dom.det(r,s)';
drdx  = 'dom.drdx(r,s)';
drdy  = 'dom.drdy(r,s)';
dsdx  = 'dom.dsdx(r,s)';
dsdy  = 'dom.dsdy(r,s)';
drdxx = 'dom.drdxx(r,s)';
drdxy = 'dom.drdxy(r,s)';
drdyy = 'dom.drdyy(r,s)';
dsdxx = 'dom.dsdxx(r,s)';
dsdxy = 'dom.dsdxy(r,s)';
dsdyy = 'dom.dsdyy(r,s)';

dxx  = op_xy.dxx;
dxy  = op_xy.dxy;
dyy  = op_xy.dyy;
dx   = op_xy.dx;
dy   = op_xy.dy;
b_xy = op_xy.b;

% d/dxx term
if ( isa( dxx, 'function_handle') )
    dxx = 'dxx(dom.x(r,s), dom.y(r,s))';
else
    dxx = num2str(dxx);
end
drr = ['@(r,s) ' dxx '.*' det '.*' drdx '.^2'];
drs = ['@(r,s) 2*' dxx '.*' det '.*' drdx '.*' dsdx ''];
dss = ['@(r,s) ' dxx '.*' det '.*' dsdx '.^2'];
dr  = ['@(r,s) ' dxx '.*' drdxx ''];
ds  = ['@(r,s) ' dxx '.*' dsdxx ''];

% d/dxy term
if ( isa( dxy, 'function_handle') )
    dxy = 'dxy(dom.x(r,s), dom.y(r,s))';
else
    dxy = num2str(dxy);
end
drr = [drr ' + ' dxy '.*' det '.*' drdx '.*' drdy ''];
drs = [drs ' + 2*' dxy '.*' det '.*(' drdy '.*' dsdx '+' drdx '.*' dsdy ')'];
dss = [dss ' + ' dxy '.*' det '.*' dsdx '.*' dsdy ''];
dr  = [dr  ' + ' dxy '.*' drdxy ''];
ds  = [ds  ' + ' dxy ' .*' dsdxy ''];

% d/dyy term
if ( isa( dyy, 'function_handle') )
    dyy = 'dyy(dom.x(r,s), dom.y(r,s))';
else
    dyy = num2str(dyy);
end
drr = [drr ' + ' dyy '.*' det '.*' drdy '.^2'];
drs = [drs ' + 2*' dyy '.*' det '.*' drdy '.*' dsdy ''];
dss = [dss ' + ' dyy '.*' det '.*' dsdy '.^2'];
dr  = [dr  ' + ' dyy '.*' drdyy ''];
ds  = [ds  ' + ' dyy '.*' dsdyy ''];

% d/dx term
if ( isa( dx, 'function_handle') )
    dx = 'dx(dom.x(r,s), dom.y(r,s))';
else
    dx = num2str(dx);
end
dr = [dr ' + ' dx '.*' det '.^2.*' drdx ''];
ds = [ds ' + ' dx '.*' det '.^2.*' dsdx ''];

% d/dy term
if ( isa( dy, 'function_handle') )
    dy = 'dy(dom.x(r,s), dom.y(r,s))';
else
    dy = num2str(dy);
end
dr = [dr ' + ' dy '.*' det '.^2.*' drdy ''];
ds = [ds ' + ' dy '.*' det '.^2.*' dsdy ''];

drr = eval(drr);
drs = eval(drs);
dss = eval(dss);
dr  = eval(dr);
ds  = eval(ds);

if ( isa(dom, 'ultraSEMTri') )
    sing = @(r,s) (r/2-1/2).^2;
    syms r s
    e = simplify( sing(r,s) .* drr(r,s) ); drr = @(r,s) eval(e);
    e = simplify( sing(r,s) .* drs(r,s) ); drs = @(r,s) eval(e);
    e = simplify( sing(r,s) .* dss(r,s) ); dss = @(r,s) eval(e);
    e = simplify( sing(r,s) .* dr(r,s) );  dr  = @(r,s) eval(e);
    e = simplify( sing(r,s) .* ds(r,s) );  ds  = @(r,s) eval(e);

    % Zeroth term
    if ( isa(b_xy, 'function_handle') )
        b_rs = @(r,s) sing(r,s) .* dom.scl(r,s) .* b_xy(dom.x(r,s), dom.y(r,s));
    else
        b_rs = @(r,s) sing(r,s) .* dom.scl(r,s) .* b_xy;
    end

    % Righthand side
    if ( isa(rhs_xy, 'function_handle') )
        rhs_rs = @(r,s) sing(r,s) .* dom.scl(r,s) .* rhs_xy(dom.x(r,s), dom.y(r,s));
    else
        rhs_rs = @(r,s) sing(r,s) .* dom.scl(r,s) .* rhs_xy;
    end
else
    % Zeroth term
    if ( isa(b_xy, 'function_handle') )
        b_rs = @(r,s) dom.scl(r,s) .* b_xy(dom.x(r,s), dom.y(r,s));
    else
        b_rs = @(r,s) dom.scl(r,s) .* b_xy;
    end

    % Righthand side
    if ( nargin > 2 )
        if ( isa(rhs_xy, 'function_handle') )
            rhs_rs = @(r,s) dom.scl(r,s) .* rhs_xy(dom.x(r,s), dom.y(r,s));
        else
            rhs_rs = @(r,s) dom.scl(r,s) .* rhs_xy;
        end
    end
end

op_rs = ultraSEMPDO({drr, drs, dss}, {dr, ds}, b_rs);

end
