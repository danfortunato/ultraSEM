function normSol = norm(sol, varargin)
%NORM   Norm of an ULTRASEM.SOL.
%   For ULTRASEM.SOL objects:
%       NORM(SOL) = sqrt(integral of abs(SOL)^2).
%       NORM(SOL, 2) is the same as NORM(SOL).
%       NORM(SOL, 'fro') is also the same as NORM(SOL).
%       NORM(SOL, 1) = integral of abs(SOL).
%       NORM(SOL, P) = (integral of abs(SOL)^P)^(1/P).
%       NORM(SOL, inf) = estimated global maximum in absolute value.
%       NORM(SOL, -inf) = estimated global minimum in absolute value.
%       NORM(SOL, 'max') is the same as NORM(SOL, inf).
%       NORM(SOL, 'min') is the same as NORM(SOL, -inf).
%       NORM(SOL, 'H1') = sqrt( ||u||_2^2 + ||grad(u)||_2^2 )
%       NORM(SOL, 'H2') = sqrt( ||u||_H1^2 + ||grad(diffx(u))||_2^2
%                                          + ||grad(diffy(u))||_2^2 )
%       NORM(SOL, 'lap') = sqrt( ||u||_2 + ||lap(u)||_2^2 )
%
%   NORM(SOL, 'all') and NORM(SOL, P, 'all') return an array of norms on
%   each patch of SOL.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

% Parse arguments.
p = 2;
reduce = true;
if ( nargin == 2 )
    if ( strcmp(varargin{1}, 'all') )
        reduce = false;
    else
        p = varargin{1};
    end
elseif ( nargin == 3 )
    p = varargin{1};
    if ( strcmp(varargin{2}, 'all') )
        reduce = false;
    end
end

% Empty ULTRASEM.SOL has norm 0.
if ( isempty(sol) )
    normSol = 0;
    return
end

switch ( p )
    case 'fro'
        normSol = norm(sol, 2);
        reduceFun = @(x) sqrt( sum(x.^2) );

    case {inf, 'inf', 'max'}
        normSol = cellfun(@(u) max(abs(u(:))), coeffs2vals(sol.coeffs));
        reduceFun = @max;

    case {-inf, '-inf', 'min'}
        normSol = cellfun(@(u) min(abs(u(:))), coeffs2vals(sol.coeffs));
        reduceFun = @min;

    case 'H1'
        [sol_x, sol_y] = grad(sol);
        norm_L2 = norm(sol,   'all');
        norm_x  = norm(sol_x, 'all');
        norm_y  = norm(sol_y, 'all');
        normSol = sqrt( norm_L2.^2 + norm_x.^2 + norm_y.^2 );
        reduceFun = @(x) sqrt( sum(x.^2) );

    case 'H2'
        [sol_x, sol_y]   = grad(sol);
        [sol_xx, sol_xy] = grad(sol_x);
        [sol_yx, sol_yy] = grad(sol_y);
        norm_H1 = norm(sol, 'H1', 'all');
        norm_xx = norm(sol_xx, 'all'); norm_xy  = norm(sol_xy, 'all');
        norm_yx = norm(sol_yx, 'all'); norm_yy  = norm(sol_yy, 'all');
        normSol = sqrt( norm_H1.^2 + norm_xx.^2 + norm_xy.^2 + ...
                                     norm_yx.^2 + norm_yy.^2 );
        reduceFun = @(x) sqrt( sum(x.^2) );

    case 'lap'
        norm_L2  = norm(sol,      'all');
        norm_lap = norm(lap(sol), 'all');
        normSol = sqrt( norm_L2.^2 + norm_lap.^2 );
        reduceFun = @(x) sqrt( sum(x.^2) );

    otherwise
        if ( isnumeric(p) && isreal(p) )
            if ( abs(round(p) - p) < eps )
                normSol = integrate(sol, p);
                reduceFun = @(x) ( sum(x.^p) ).^(1/p);
            else
                error('ULTRASEM:SOL:norm:norm', ...
                    'ULTRASEM.SOL does not support this norm.');
            end
        else
            error('ULTRASEM:SOL:norm:unknown', 'Unknown norm.');
        end

end

% Combine norms on each element.
if ( reduce )
    normSol = reduceFun(normSol);
end

end

function int = integrate(sol, p)
%INTEGRATE   Compute (integral of abs(SOL)^P)^(1/P)

% P should be an integer.
p = round(p);

% If a patch uses an N x N discretization, then quadrature is performed on
% that patch using N*P points.
int = zeros(length(sol),1);
for k = 1:length(sol)
    [ny,nx] = size(sol.coeffs{k});
    qx = nx*p; qy = ny*p;
    U = zeros(qy,qx);
    U(1:ny,1:nx) = sol.coeffs{k};
    V = util.coeffs2vals( util.coeffs2vals(U).' ).';
    wx = util.quadwts(qx); wx = wx(:);
    wy = util.quadwts(qy); wy = wy(:);
    int(k) = sum(sum(abs(V).^p .* wy .* wx.')).^(1/p);
end

end
