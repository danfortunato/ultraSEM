function P = initialize(dom, op, rhs, n)
%INITIALIZE   Initialize an array of ultraSEMLeaf objects.
%
%   P = ULTRASEMLEAF.INITIALIZE(DOM, OP, RHS) returns a cell array P of
%   ULTRASEMLEAF objects which contain the solution and D2N operators for
%   the PDO specified by OP (which may be an ultraSEMPDO or a cell array)
%   on the domain DOM (which may be an ultraSEMDomain or an Mx4 matrix
%   containing the coordinates of each patch) with zero righthand side.
%
%   P = ULTRASEMLEAF.INITIALIZE(DOM, OP, RHS, N) specifies the
%   discretization size N on each patch. When not given (or empty) N
%   defaults to 21.

% Copyright 2018 by Nick Hale and Dan Fortunato.

% Default discretization size:
default_n = 21;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%% PARSE INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ( nargin < 3 )
    % Default to homogeneous RHS:
    rhs = 0;
end
if ( nargin < 4 || isempty(n) )
    % Default value of n:
    n = default_n;
end
if ( isa(dom, 'ultraSEMDomain') )
    % Typically we are given a tree structure. Convert to array:
    dom = dom.domain;
end
numPatches = size(dom, 1);
if ( ~isa(op, 'ultraSEMPDO') )
    % PDE given as cell. Convert to PDO:
    op = ultraSEMPDO(op);
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check for a constant-coefficient PDO on a uniform domain:
[~, isConstant] = feval(op, dom(1,1), dom(1,3));
if ( ~isConstant || any(diff(diff(dom(:,1:2),1,2))) || ...
                    any(diff(diff(dom(:,3:4),1,2))) )
    error('ULTRASEM:ULTRASEMLEAF:initialize:nonconst', ...
        'Cannot handle variable-coefficient PDEs or non-rectangular domains yet.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%% DEFINE CORNER CONVENTIONS %%%%%%%%%%%%%%%%%%%%%%%%%
numIntDOF = (n-2)^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%% DEFINE BOUNDARY NODES %%%%%%%%%%%%%%%%%%%%%%%%%%%%
[XX, YY] = chebpts2(n); % Chebyshev points and grid.
leftIdx  = sub2ind([n n], (1:n).', ones(n,1));
rightIdx = sub2ind([n n], (1:n).', n*ones(n,1));
downIdx  = sub2ind([n n], ones(n,1), (1:n).');
upIdx    = sub2ind([n n], n*ones(n,1), (1:n).');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%% DEFINE BOUNDARY CONDITIONS %%%%%%%%%%%%%%%%%%%%%%%%%

% Construct normal derivatives conditions along the four edges:
I = speye(n);
lbc_d = kron( (-1).^(0:n-1).*(0:n-1).^2, I );
rbc_d = kron( ones(1,n).*(0:n-1).^2, I );
dbc_d = kron( I, (-1).^(0:n-1).*(0:n-1).^2 );
ubc_d = kron( I, ones(1,n).*(0:n-1).^2 );
bcrows_d = [ lbc_d ; rbc_d ; dbc_d ; ubc_d ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%% CONSTANT RHS? %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define scalar RHSs:
if ( isnumeric(rhs) )
    if ( isscalar(rhs) )
        rhs_eval = [ rhs; zeros(numIntDOF-1,1) ];
        constantRHS = true;
    else
        % The user gave us coeffs:
        assert(all(size(rhs) == [numIntDOF numPatches]), 'Invalid RHS.');
        rhs_eval = rhs;
        constantRHS = false;
    end
else
    % We need to evaluate the RHS:
    rhs_eval = zeros(numIntDOF, numPatches);
    constantRHS = false;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%% SOLVE LOCAL PROBLEMS %%%%%%%%%%%%%%%%%%%%%%%%%

% Solution operator:
% This involves solving a O(p^2) x O(p^2) almost-banded-block-banded system
% O(p) times. Currently this is O(p^4*p) = O(p^5), but can be done in
% O(p^4) using e.g. Woodbury.
[S, A] = buildSolOp(op, dom(1,:), rhs_eval(:,1), n);

% Dirichlet-to-Neumann map:
D2N = bcrows_d * S;

% Initialize
P = cell(numPatches, 1);

% Scaling (for all patches):
domx = dom(1,1:2);
domy = dom(1,3:4);
sclx = 2/diff(domx);
scly = 2/diff(domy);

% Loop over each patch:
for k = 1:numPatches

    % Define the boundary nodes for this patch:
    x = (XX+1)/sclx + dom(k,1);
    y = (YY+1)/scly + dom(k,3);
    % Store the four sides separately:
    xy = {[ x(leftIdx)  y(leftIdx)  ] ;
          [ x(rightIdx) y(rightIdx) ] ;
          [ x(downIdx)  y(downIdx)  ] ;
          [ x(upIdx)    y(upIdx)    ] };

    % Evaluate non-constant RHSs:
    if ( ~isnumeric(rhs) )
        vals = feval(rhs, x, y);
        % Convert to coeffs:
        coeffs = chebtech2.vals2coeffs(chebtech2.vals2coeffs(vals).').';
        % Map the RHS to the right ultraspherical space.
        lmap = ultraS.convertmat(n, 0, 1);
        rmap = ultraS.convertmat(n, 0, 1);
        coeffs = lmap * coeffs * rmap.';
        coeffs = coeffs(1:n-2, 1:n-2);
        rhs_eval(:,k) = reshape(coeffs, numIntDOF, 1);
    end

    % Assemble the patch:
    P{k} = ultraSEMLeaf(dom(k,:), S, D2N, xy, n, op);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Append particular parts:
if ( ~constantRHS )
    Si = A \ rhs_eval;
    Gx = zeros(2, n); Gy = zeros(2, n);
    Bx = [(-1).^(0:n-1); ones(1,n)];
    By = [(-1).^(0:n-1); ones(1,n)];
    [By, Gy, Py] = canonicalBC(By, Gy);
    [Bx, Gx, Px] = canonicalBC(Bx, Gx);
    for k = 1:numPatches
        S = imposeBCs(Si(:,k), Px, Py, Bx, By, Gx, Gy, n);
        P{k}.S(:,end) = S;
        P{k}.D2N(:,end) = P{k}.D2N(:,end) + bcrows_d * S;
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%% DEFINE OPERATORS %%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [S, A] = buildSolOp(pdo, dom, rhs, n)
%BUILDSOLOP  Build the solution operator on a patch using chebop2 ideas.
%
%   S = solution operator on patch
%   A = discretized operator on patch

    % Convert variable coefficients of PDO to chebfun2 in order to get a
    % separable representation.
    for field = fieldnames(pdo)'
        coeff = pdo.(field{1});
        if ( isa(coeff, 'function_handle') )
            pdo.(field{1}) = chebfun2(coeff); % Domain?
        end
    end

    % Discretize the separable ODOs that make the PDO
    CC = discretizeODOs(pdo, dom, n);

    % Encode all possible BCs
    % Note: this could be precomputed for a given n.
    Bx = [(-1).^(0:n-1); ones(1,n)];
    By = [(-1).^(0:n-1); ones(1,n)];
    Gx = zeros(2, n, 4*n);
    Gy = zeros(2, n, 4*n);
    P = compatibleProjection(n);
    for k = 1:4*n
        Gx(1,:,k) = P(1:n,k);
        Gx(2,:,k) = P(n+1:2*n,k);
        Gy(1,:,k) = P(2*n+1:3*n,k);
        Gy(2,:,k) = P(3*n+1:end,k);
    end

    % Canonicalize:
    [By, Gy, Py] = canonicalBC(By, Gy);
    [Bx, Gx, Px] = canonicalBC(Bx, Gx);

    % Remove 4n-4 degrees of freedom by enforcing the boundary constraints
    [CC, BC] = eliminateBCs(CC, Bx, By, Gx, Gy, n);

    % Form the n^2 x n^2 linear system
    A = spalloc((n-2)^2, (n-2)^2, n*(n-2)^2 + (n-2)^2);
    for k = 1:size(CC, 1)
        A = A + kron(CC{k,2}, CC{k,1});
    end
    BC = reshape(BC, (n-2)^2, 4*n);

    % Solve for every possible BC
%     S22 = A \ [BC, rhs];
    S22 = schurSolve(A, [BC, rhs], 2*n-2);

    % Add in the boundary data
    Gx(:,:,end+1) = zeros(2, n);
    Gy(:,:,end+1) = zeros(2, n);
    S = imposeBCs(S22, Px, Py, Bx, By, Gx, Gy, n);

end

function CC = discretizeODOs(pdo, dom, n)
%DISCRETIZEODOS   Discretize the PDO into a separable representation.

    terms = fieldnames(pdo);
    CC = cell(size(terms,1), 1);

    for k = 1:size(terms,1)

        % Get the PDO coefficient and the differentiation order in x and y.
        coeff = pdo.(terms{k});
        dx = length(strfind(terms{k}, 'x'));
        dy = length(strfind(terms{k}, 'y'));

        % Check if the coefficient is zero
        if ( (isscalar(coeff) && coeff ~= 0) || isa(coeff, 'chebfun2') )

            % Define operators
            CC{k} = cell(length(coeff), 2);
            Sx = ultraS.convertmat(n, dx, 1);
            Sy = ultraS.convertmat(n, dy, 1);
            Dx = (2/diff(dom(3:4)))^dx * ultraS.diffmat(n, dx);
            Dy = (2/diff(dom(1:2)))^dy * ultraS.diffmat(n, dy);

            if ( isscalar(coeff) )
                % Constant coefficient
                CC{k}{1,1} = coeff .* Sx * Dx;
                CC{k}{1,2} =          Sy * Dy;
            else
                % Variable coefficient
                [C, D, R] = cdr(coeff);
                % Make a multiplication operator for each slice of chebfun2
                for r = 1:length(coeff)
                    Mx = ultraS.multmat( n, D(r,r) * R(:,r), dx );
                    My = ultraS.multmat( n,          C(:,r), dy );
                    CC{k}{r,1} = Sx * Mx * Dx;
                    CC{k}{r,2} = Sy * My * Dy;
                end
            end
        end
    end

    % Flatten into a rank x 2 cell array
    CC = cat(1, CC{:});

end

function [CC, BC] = eliminateBCs(CC, Bx, By, Gx, Gy, n)
%ELIMINATEBCS   Remove degrees of freedom from solution due to BCs.

    BC = zeros(n, n, 4*n);
    for k = 1:size(CC, 1)
        [CC{k,1}, BC] = zeroDOF(CC{k,1}, CC{k,2}, BC, By, Gy);
        [CC{k,2}, BC] = zeroDOF(CC{k,2}, CC{k,1}, BC, Bx, Gx, true);
    end

    % Remove degrees of freedom.
    for k = 1:size(CC, 1)
        CC{k,1} = CC{k,1}(1:n-2, 3:n);
        CC{k,2} = CC{k,2}(1:n-2, 3:n);
    end
    % Truncation of righthand side.
    BC = BC(1:n-2, 1:n-2, :);

end

function [C1, E] = zeroDOF(C1, C2, E, B, G, trans)
%ZERODOF   Eliminate so degrees of freedom in the matrix equation can be
%removed.

    if ( nargin < 6)
        trans = false;
    end

    for ii = 1:size(B, 1) % For each boundary condition, zero a column.
        for kk = 1:size(C1, 1)
            if ( abs(C1(kk,ii)) > 10*eps )
                c = C1(kk,ii); % Constant required to zero entry out.
                C1(kk,:) = C1(kk,:) - c*B(ii,:);
                for ll = 1:size(G, 3) % Do this for every BC.
                    if ( trans )
                        E(:,kk,ll) = E(:,kk,ll).' - c*G(ii,:,ll)*C2.';
                    else
                        E(kk,:,ll) = E(kk,:,ll) - c*G(ii,:,ll)*C2.';
                    end
                end
            end
        end
    end

end

function [B, G, P] = canonicalBC(B, G)
%CANONICALBC   Form a linear combination of the boundary conditions
%so that they can be used for imposing on the PDE.

    P = nonsingularPermute(B);
    B = B*P;
    [L, B] = lu(B);

    % Scale so that B is unit upper triangular.
    if ( min(size(B)) > 1 )
        D = diag(1./diag(B));
    elseif ( ~isempty(B) )
        D = 1./B(1,1);
    else
        D = []; % No boundary conditions.
    end
    B = D*B;

    for k = 1:size(G,3)
        G(:,:,k) = L \ G(:,:,k);
        G(:,:,k) = D*G(:,:,k);
    end

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
        error('ULTRASEM:buildSolOp:nonsingularPermute:BCs', ...
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

function P = compatibleProjection(n)
%COMPATIBLEPROJECTION   Construct an orthogonal projection onto compatible
%boundary conditions.
%
%   P = 4n x 4n matrix that imposes continuity at the corners.

    % Order is: [lbc; rbc; dbc; ubc]

    % Top left corner
    A = zeros(4,4*n);
    A(1,1:n) = ones(1,n);
    A(1,3*n+1:end) = -(-1).^(0:n-1);

    % Bottom left corner
    A(2,1:n) = (-1).^(0:n-1);
    A(2,2*n+1:3*n) = -(-1).^(0:n-1);

    % Top right corner
    A(3,n+1:2*n) = ones(1,n);
    A(3,3*n+1:end) = -ones(1,n);

    % Bottom right corner
    A(4,n+1:2*n) = (-1).^(0:n-1);
    A(4,2*n+1:3*n) = -ones(1,n);

    % Build the projection:
    VV = null(A);
    P = VV * VV';

end

function x = schurSolve(A, b, m)
% Fast solution of A*x = b where A is banded + m dense rows via Schur
% complement factorisation.

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
