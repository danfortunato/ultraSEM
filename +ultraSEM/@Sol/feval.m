function coeffs = feval(sol, x, y)
%FEVAL   Evaluate an ULTRASEM.SOL at one or more points.
%   FEVAL(SOL, X, Y) evaluates SOL at the point(s) (X, Y).

if ( ~isnumeric(sol.domain) )
    coeffs = fevalMapped(sol, x, y);
    return
end

minx = min(sol.domain(:,1));
miny = min(sol.domain(:,3));

sx = size(x);
coeffs = 0*x(:);
x = x(:);
y = y(:);

for k = 1:length(sol)

    dom = sol.domain(k,:);
    idx = (x > dom(1) | (dom(1) == minx & x == minx)) & x <= dom(2) & ...
          (y > dom(3) | (dom(3) == miny & y == miny)) & y <= dom(4);

    if ( any(idx) )
        % Map x and y to [-1,1]^2:
        sclx = 2/diff(dom(1:2));
        scly = 2/diff(dom(3:4));
        xm = sclx*(x(idx)-dom(1))-1;
        ym = scly*(y(idx)-dom(3))-1;
        coeffs(idx) = util.clenshaw2d(sol.coeffs{k}, xm, ym);
    end

end

coeffs = reshape(coeffs, sx);

end

function u = fevalMapped(sol, x, y)
%FEVALMAPPED   Evaluate an ULTRASEM.SOL at one or more points.
%   FEVALMAPPED(SOL, X, Y) evaluates SOL at the point(s) (X, Y), where the
%   domain of SOL consists of a generic ULTRASEM.MAPPING.

sx = size(x);
u = 0*x(:);
x0 = x(:);
y0 = y(:);

for k = 1:size(sol.domain, 1)

    map = sol.domain(k,:);
    % Convert to single to avoid mapping issues
    x = map.r(x0, y0);
    y = map.s(x0, y0);
    dom = [-1 1 -1 1];

    xs = single(x);
    ys = single(y);

    idx = ( xs == real(xs) & ys == real(ys) & ...
            xs >= dom(1) & xs <= dom(2)     & ...
            ys >= dom(3) & ys <= dom(4) );

    if ( any(idx) )
        u(idx) = util.clenshaw2d(sol.coeffs{k}, x(idx), y(idx));
    end

end

u = reshape(u, sx);

end
