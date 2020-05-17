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
%
%   NORM(SOL, 'all') and NORM(SOL, P, 'all') return a cell array of norms
%   on each patch of SOL.

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
        reduceFun = @sum;

    case {inf, 'inf', 'max'}
        normSol = cellfun(@(u) max(abs(u(:))), coeffs2vals(sol.u));
        normSol = num2cell(normSol);
        reduceFun = @max;

    case {-inf, '-inf', 'min'}
        normSol = cellfun(@(u) min(abs(u(:))), coeffs2vals(sol.u));
        normSol = num2cell(normSol);
        reduceFun = @min;

    otherwise
        if ( isnumeric(p) && isreal(p) )
            if ( abs(round(p) - p) < eps )
                normSol = integrate(sol, p);
                reduceFun = @sum;
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
    normSol = reduceFun(cell2mat(normSol));
end

end

function int = integrate(sol, p)
%INTEGRATE   Compute (integral of abs(SOL)^P)^(1/P)

% P should be an integer.
p = round(p);

% If a patch uses an N x N discretization, then quadrature is performed on
% that patch using N*P points.
int = cell(length(sol),1);
for k = 1:length(sol)
    [ny,nx] = size(sol.u{k});
    qx = nx*p; qy = ny*p;
    U = zeros(qy,qx);
    U(1:ny,1:nx) = sol.u{k};
    V = util.coeffs2vals( util.coeffs2vals(U).' ).';
    wx = util.quadwts(qx); wx = wx(:);
    wy = util.quadwts(qy); wy = wy(:);
    int{k} = sum(sum(abs(V).^p .* wy .* wx.')).^(1/p);
end

end
