function P = updateRHS(P, rhs)
%UPDATERHS   Update RHS of an ULTRASEM.LEAF object.
%   P = UPDATERHS(P, F) replaces the existing RHS of an initialized
%   ULTRASEM.LEAF object P with that given in F, which may be a constant or
%   a function handle. Since only one subproblem needs to be solved on each
%   patch in the update process (rather than O(n) in the original
%   initialization) this can lead to considerable performance gains when
%   solving for multiple RHSs.
%
% See also ULTRASEM.LEAF.INITIALIZE.

% Copyright 2018 by Nick Hale and Dan Fortunato.

% Developer note: At user levels this is typically called with 
%  >> P.rhs = F

p = P.p;
dom = P.domain;
numIntDOF = (p-2)^2;

% TODO: This is unfortuate...
op = ultraSEM.PDO({0,0,0}); % Entries don't matter.
[~, rhs] = transformPDO(dom, op, rhs);

% Define scalar RHSs:
if ( isnumeric(rhs) )
    if ( isscalar(rhs) )
        % Constant RHS.
        rhs = [ rhs; zeros(numIntDOF-1,1) ];
    end
else
    % Evaluate non-constant RHS:
    [X, Y] = util.chebpts2(p); % Chebyshev points and grid.
    rhs = evaluateRHS(rhs, X, Y, p, numIntDOF);
end

Ainv = P.Ainv;
if ( isempty(Ainv) )
    error('ULTRASEM:leaf:updateRHS:operatorNotStored', ...
        'Discretized operator A was not stored. Cannot update RHS.');
    % TODO: Perhaps we can store A _OR_ the PDO. In the latter case, rebuild A.
end
Sp = Ainv(rhs);

% Encode BCs:
Gx = zeros(2, p); Gy = zeros(2, p);
Bx = [(-1).^(0:p-1); ones(1,p)];
By = [(-1).^(0:p-1); ones(1,p)];
[By, Gy, Py] = canonicalBC(By, Gy);
[Bx, Gx, Px] = canonicalBC(Bx, Gx);
Sp = imposeBCs(Sp, Px, Py, Bx, By, Gx, Gy, p);
P.S(:,end) = Sp;

% Normal derivative:
normal_d = transformNormalD(dom, p);
P.D2N(:,end) = normal_d * Sp;   

end

function out = evaluateRHS(rhs, x, y, p, numIntDOF)
vals = feval(rhs, x, y);
% Convert to coeffs:
coeffs = util.vals2coeffs(util.vals2coeffs(vals).').';
% Map the RHS to the right ultraspherical space:
lmap = util.convertmat(p, 0, 1);
rmap = util.convertmat(p, 0, 1);
coeffs = lmap * coeffs * rmap.';
coeffs = coeffs(1:p-2, 1:p-2);
out = reshape(coeffs, numIntDOF, 1);
end

function [B, G, P] = canonicalBC(B, G)
%CANONICALBC   Form a linear combintation of the boundary conditions
%so that they can be used for imposing on the PDE.

P = nonsingularPermute(B);
B = B*P;
[L, B] = lu(B);
G = L \ G;

% Scale so that B is unit upper triangular.
if ( min(size(B)) > 1 )
    D = diag(1./diag(B));
elseif ( ~isempty(B) )
    D = 1./B(1,1);
else
    D = []; % No boundary conditions.
end
B = D*B;
G = D*G;

end

function P = nonsingularPermute(B)
%NONSINGULARPERMUTE   Permute the columns of B to ensure that the principal
%m*m submatrix of B is nonsingular, where m = size(B, 1).
%
% Note: This is needed for solving the matrix equations with linear
% constraints, see DPhil thesis of Alex Townsend (section 6.5).

m = size(B, 1);
k = 1;

% [TODO]: improve this check.
% Try each mxm block in a linear fashion:
while ( rank(B(:,k:m+k-1)) < m )
    k = k+1;
    if ( m+k > size(B, 2) )
        error('CHEBFUN:CHEBOP2:discretize:nonsingularPermute:BCs', ...
            'Boundary conditions are linearly dependent.');
    end
end

P = speye(size(B, 2));
P = P(:,[k:m+k-1, 1:k-1, m+k:end]);

end

function S = imposeBCs(S22, Px, Py, Bx, By, Gx, Gy, n)
%IMPOSEBCS   Impose the boundary conditions on the solution.

    By = By * Py;
    Bx = Bx * Px;

    nc = size(Gx, 3);
    S = zeros(n^2, nc);
    for k = 1:nc
        % Recombine the boundary conditions.
        X22 = reshape(S22(:,k), n-2, n-2);
        X12 = By(:,1:2) \ (Gy(:,3:n,k) - By(:,3:n)*X22);
        X = [ X12; X22 ];
        X2 = Bx(:,1:2) \ (Gx(:,1:n,k) - Bx(:,3:n)*X.');
        X = [ X2.' X ];
        S(:,k) = X(:);
    end

end