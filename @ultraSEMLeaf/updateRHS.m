function P = updateRHS(P, rhs)
%UPDATERHS   Update RHS of an ultraSEMLeaf object.
%   P = UPDATERHS(P, F) replaces the existing RHS of an initialized
%   ultraSEMLeaf object P with that given in F, which may be a constant or
%   a function handle. Since only one subproblem needs to be solved on each
%   patch in the update process (rather than O(n) in the original
%   initialization) this can lead to considerable performance gains when
%   solving for multiple RHSs.
%
% See also ULTRASEMLEAF.INITIALIZE.

% Copyright 2018 by Nick Hale and Dan Fortunato.

Ainv = P.Ainv;
if ( isempty(Ainv) )
    error('ULTRASEM:ULTRASEMLEAF:updateRHS:operatorNotStored', ...
        'Discretised operator A was not stored. Cannot update RHS.');
    % TODO: Perhaps we can store A _OR_ the PDO. In the latter case, rebuild A.
end
n = P.n;
dom = P.domain;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%% DEFINE GRID ON [-1 1] AND SET INDICIES %%%%%%%%%%%%%%%%%%%%%%%%%%%
numIntDOF = (n-2)^2;
[XX, YY] = util.chebpts2(n);

mydom = dom;
if ( isRect(dom) )
    %TODO: Fix this hack.
    dom = quad2rect(dom.v);
end

% Define the scaling for this domain:
if ( isnumeric(dom) )
    domx = dom(1:2);
    domy = dom(3:4);
else
    domx = dom.domain(1:2);
    domy = dom.domain(3:4);
end
% Define the grid for this patch:
sclx = 2/diff(domx);
scly = 2/diff(domy);
x = (XX+1)/sclx + domx(1);
y = (YY+1)/scly + domy(1);

% Evaluate non-constant RHSs:
if ( isnumeric(rhs) && isscalar(rhs) )
    rhs = [ rhs; zeros(numIntDOF-1,1) ];
elseif ( ~isnumeric(rhs) )
    vals = feval(rhs, x, y);
    % Convert to coeffs:
    coeffs = util.vals2coeffs(util.vals2coeffs(vals).').';
    % Map the RHS to the right ultraspherical space:
    lmap = util.convertmat(n, 0, 1);
    rmap = util.convertmat(n, 0, 1);
    coeffs = lmap * coeffs * rmap.';
    coeffs = coeffs(1:n-2, 1:n-2);
    rhs = reshape(coeffs, numIntDOF, 1);
end

% Solve with the new RHS:
% S = schurSolve(A, rhs, 2*n-2);
S = Ainv(rhs);

% Impose zero Dirichlet BCs:
Gx = zeros(2, n); Gy = zeros(2, n);
Bx = [(-1).^(0:n-1); ones(1,n)];
By = [(-1).^(0:n-1); ones(1,n)];
[By, Gy, Py] = canonicalBC(By, Gy);
[Bx, Gx, Px] = canonicalBC(Bx, Gx);
S = imposeBCs(S, Px, Py, Bx, By, Gx, Gy, n);

% Amend final column of the solution operator:
P.S(:,end) = S;

% Normal derivative operator:
if ( ~isRect(mydom) )
    normal_d = transformNormalD(mydom, P.xy, n);
else
    % Construct normal derivatives conditions along the four edges:
    I = speye(n);
    lbc_d = sclx*kron( (-1).^(0:n-1).*(0:n-1).^2, I );
    rbc_d = sclx*kron( ones(1,n).*(0:n-1).^2, I );
    dbc_d = scly*kron( I, (-1).^(0:n-1).*(0:n-1).^2 );
    ubc_d = scly*kron( I, ones(1,n).*(0:n-1).^2 );
    normal_d = [ lbc_d ; rbc_d ; dbc_d ; ubc_d ];
end

% Amend final column of the D2N operator:
P.D2N(:,end) = normal_d * S;

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

function x = schurSolve(A, b, m)
%SCHURSOLVE   Fast solution of A*x = b where A is banded + m dense rows via
%Schur complement factorisation.

    doRowScaling = true;

    na = size(A,2);
    nb = size(b,2);
    if ( na <= m )
        x = A\b;
        return
    end

    i1 = 1:m;
    i2 = m+1:na;
    i3 = nb+(1:m);

    if ( doRowScaling )
        % Row scaling to improve accuracy
        AA = A(i2,i2);
        s = 1./ max(1, max(abs(AA), [], 2) );
        AA = bsxfun(@times, s, AA);
        bb = s.*[b(i2,:), A(i2,i1)];
        c = AA\bb;
    else
        c = A(i2,i2)\[b(i2,:), A(i2,i1)];
    end

    x = (A(i1,i1) - A(i1,i2)*c(:,i3)) \ (b(i1,:) - A(i1,i2)*c(:,1:nb));
    y = c(:,1:nb) - c(:,i3)*x;
    x = [x ; y];

end
