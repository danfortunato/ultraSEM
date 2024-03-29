function L = initialize(dom, op, varargin)
%INITIALIZE   Initialize an array of ULTRASEM.LEAF objects.
%   L = ULTRASEM.LEAF.INITIALIZE(DOM, OP) returns a cell array L of
%   ULTRASEM.LEAF objects which contain the solution and D2N operators for
%   the PDO specified by OP (which may be an ULTRASEM.PDO or a cell array)
%   on the domain DOM (which may be an ULTRASEM.MAPPING) with zero
%   righthand side.
%
%   L = ULTRASEM.LEAF.INITIALIZE(DOM, OP, RHS) is as above, but with the
%   righthand side RHS, which may be a scalar or a function handle.
%
%   L = ULTRASEM.LEAF.INITIALIZE(DOM, OP, RHS, P) specifies the
%   discretization size P on each patch. When not given (or empty), P
%   defaults to ULTRASEM.PREF.DISCSIZE. P may be a scalar (in which case
%   the same discretization size is used on each patch), or a vector of
%   length LENGTH(DOM) (in which case a different disctization size is used
%   on each patch.
%
%   L = ULTRASEM.LEAF.INITIALIZE(..., PREF) uses the preferences specified
%   in the ULTRASEM.PREF object PREF. (See ULTRASEM.PREF for details on the
%   various preference options and their defaults.)

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%% PARSE INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set defaults for unspecified arguments:
rhs = 0; % Default to homogeneous problem.
p = [];
pref = ultraSEM.Pref();

if ( nargin == 3 )
    if ( isa(varargin{1}, 'ultraSEM.Pref') )
        pref = varargin{1}; % INITIALIZE(DOM, OP, PREF)
    else
        rhs = varargin{1};  % INITIALIZE(DOM, OP, RHS)
    end
elseif ( nargin == 4 )
    rhs = varargin{1};
    if ( isa(varargin{2}, 'ultraSEM.Pref') )
        pref = varargin{2}; % INITIALIZE(DOM, OP, RHS, PREF)
    else
        p = varargin{2};    % INITIALIZE(DOM, OP, RHS, N)
    end
elseif ( nargin == 5 )
    rhs = varargin{1};      % INITIALIZE(DOM, OP, RHS, N, PREF)
    p = varargin{2};
    pref = varargin{3};
end

assert(isa(dom, 'ultraSEM.Mapping'), 'Invalid DOM');
numPatches = size(dom, 1);

% Determine p:
if ( isempty(p) )
    % Default discretization size:
    p = pref.discSize;
elseif ( (isvector(p) && ~isscalar(p)) || isa(p, 'function_handle') )
    % Elements have varying p:
    L = cell(numPatches, 1);
    if ( isa(p, 'function_handle') )
        for j = 1:numPatches
            c = centroid(dom(j,:));
            pj = p(c(1),c(2));
            pj = max(floor(pj), 3);
            L(j) = ultraSEM.Leaf.initialize(dom(j,:), op, rhs, pj, pref);
        end
    else
        assert(numel(p) == numPatches, ...
            'Number of p''s must equal number of patches.');
        for j = 1:numPatches
            L(j) = ultraSEM.Leaf.initialize(dom(j,:), op, rhs, p(j), pref);
        end
    end
    return
else
    assert(isscalar(p), 'Invalid p.');
end

if ( ~isa(op, 'ultraSEM.PDO') )
    % PDE given as cell. Convert to PDO:
    op = ultraSEM.PDO(op);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%% DEFINE CORNER CONVENTIONS %%%%%%%%%%%%%%%%%%%%%%%%%
numIntDOF = (p-2)^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%% DEFINE REFERENCE GRID %%%%%%%%%%%%%%%%%%%%%%%%%%%%
[X, Y] = util.chebpts2(p); % Chebyshev points and grid.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%% CONSTANT RHS? %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define scalar RHSs:
if ( isnumeric(rhs) )
    if ( isscalar(rhs) )
        % Constant RHS.
        rhs_eval = zeros(numIntDOF, numPatches);
        rhs_eval(1,:) = rhs;
        constantRHS = true;
    else
        % The user gave us coeffs. (This situation is probably rare.)
        assert(all(size(rhs) == [numIntDOF numPatches]), 'Invalid RHS.');
        rhs_eval = rhs;
        constantRHS = false;
    end
else
    % We need to evaluate the RHS. Initialise to zero for now:
    rhs_eval = zeros(numIntDOF, numPatches);
    constantRHS = false;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%% SOLVE LOCAL PROBLEMS %%%%%%%%%%%%%%%%%%%%%%%%%

% Check for a constant-coefficient PDO on a uniform domain:
if ( isUniformOp(op, dom) )
    % Exploit that each patch has similar solution operator:
    L = buildConstOp();
else
    % Each patch needs a different solution operator:
    L = buildNonconstOp();
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function L = buildConstOp()
        % In the simple case where the coefficients of the PDO are constant
        % (after transformation) the solution operators on each patch are
        % (almost) the same and we can save a lot of time.
        
        % Initialize
        L = cell(numPatches, 1);
        
        % Solution operator (same for each patch):
        [S, Ainv] = buildSolOp(op, dom(1), rhs_eval(:,1), p, pref);
        
        % Dirichlet-to-Neumann map (same for each patch):
        normal_d = transformNormalD(dom(1), p);
        D2N = normal_d * S;
        
        % Loop over patches. Transform grid and evaluate RHS where req'd.
        for k = 1:numPatches
            
            % Define the domain edges for this patch:
            domk = dom(k);
            edges = [domk.v([1,2,1,4],:), domk.v([4,3,2,3],:), repmat(p,4,1)];
            
            % Evaluate non-constant RHSs:
            if ( ~isnumeric(rhs) )
                [x, y] = transformGrid(domk, X, Y);
                rhs_eval(:,k) = evaluateRHS(rhs, x, y, p, numIntDOF);
            end
            
            % Assemble the patch:
            L{k} = ultraSEM.Leaf(domk, S, D2N, edges, Ainv, normal_d);
        end
        
        % Append particular parts if necessary (i.e., non constant RHS):
        if ( ~constantRHS )
            Sp = Ainv(rhs_eval);
            % Encode BCs:
            Gx = zeros(2, p); Gy = zeros(2, p);
            Bx = [(-1).^(0:p-1); ones(1,p)];
            By = [(-1).^(0:p-1); ones(1,p)];
            [By, Gy, Py] = canonicalBC(By, Gy);
            [Bx, Gx, Px] = canonicalBC(Bx, Gx);
            for k = 1:numPatches
                Spk = imposeBCs(Sp(:,k), Px, Py, Bx, By, Gx, Gy, p);
                L{k}.S(:,end) = Spk;
                % TODO: Which of these is correct?
%                 L{k}.D2N(:,end) = L{k}.D2N(:,end) + normal_d * Spk;
                L{k}.D2N(:,end) = normal_d * Spk;
            end
        end
        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function L = buildNonconstOp()
        % In the generic case, we must fully construct the solution
        % operator on each patch.
        
        % Initialize
        L = cell(numPatches, 1);
        
        % Loop over each patch:
        for k = 1:numPatches
            
            % Define the domain edges for this patch:
            domk = dom(k,:);
            edges = [domk.v([1,2,1,4],:), domk.v([4,3,2,3],:), repmat(p,4,1)];
            
            % Transform the equation:
            [op_k, rhs_k] = transformPDO(domk, op, rhs);

            % Transform the normal derivative:
            normal_d = transformNormalD(domk, p);
            
            % Evaluate non-constant RHSs if required:
            if ( ~isnumeric(rhs_k) )
                rhs_eval(:,k) = evaluateRHS(rhs_k, X, Y, p, numIntDOF);
            end
            
            % Solution operator:
            [S, Ainv] = buildSolOp(op_k, domk, rhs_eval(:,k), p, pref);

            
            % Dirichlet-to-Neumann map:
            D2N = normal_d * S;
            
            % Assemble the patch:
            L{k} = ultraSEM.Leaf(domk, S, D2N, edges, Ainv, normal_d);
            
        end
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%% DEFINE OPERATORS %%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [S, Ainv] = buildSolOp(pdo, dom, rhs, p, pref)
%BUILDSOLOP  Build the solution operator on a patch using chebop2 ideas.
%
%   S = solution operator on patch
%   A = discretized operator on patch

% TODO: This breaks encapsulation.
if ( isRect(dom) )
    dom = rectVertices(dom);
else
    dom = [-1 1 -1 1];
end

% Discretize the separable ODOs that make the PDO
CC = discretizeODOs(pdo, dom, p);

% Encode all possible BCs
[Bx, Gx, Px, By, Gy, Py] = encodeBCs(p);

% Remove 4p-4 degrees of freedom by enforcing the boundary constraints
[CC, BC] = eliminateBCs(CC, Bx, By, Gx, Gy, p);
BC = reshape(BC, (p-2)^2, 4*p);

% Form the p^2 x p^2 linear system
A = formSystem(CC, p);

% Solve for every possible BC.
[S22, Ainv] = mysolve(A, BC, rhs, p, pref);

% Add in the boundary data
Gx(:,:,end+1) = zeros(2, p);
Gy(:,:,end+1) = zeros(2, p);
S = imposeBCs(S22, Px, Py, Bx, By, Gx, Gy, p);

end

function A = formSystem(CC, p)
A = spalloc((p-2)^2, (p-2)^2, p*(p-2)^2 + (p-2)^2);
for k = 1:size(CC, 1)
    Ak = kron(CC{k,2}, CC{k,1});
    A = A + Ak;
end
end

function [Bx, Gx, Px, By, Gy, Py] = encodeBCs(p)

persistent storage

if ( isempty(storage) || storage{1} ~= p )
    % Note: this can be precomputed for a given p.
    Bx = [(-1).^(0:p-1); ones(1,p)];
    By = [(-1).^(0:p-1); ones(1,p)];
    Gx = zeros(2, p, 4*p);
    Gy = zeros(2, p, 4*p);
    P = compatibleProjection(p);
    for k = 1:4*p
        Gx(1,:,k) = P(1:p,k);
        Gx(2,:,k) = P(p+1:2*p,k);
        Gy(1,:,k) = P(2*p+1:3*p,k);
        Gy(2,:,k) = P(3*p+1:end,k);
    end
    
    % Canonicalize:
    [By, Gy, Py] = canonicalBC(By, Gy);
    [Bx, Gx, Px] = canonicalBC(Bx, Gx);
    
    % Store for reuse.
    storage = {p, Bx, Gx, Px, By, Gy, Py};
else
    % Values are already in storage
    [~, Bx, Gx, Px, By, Gy, Py] = deal(storage{:});
    return
end

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
    if ~( (isnumeric(coeff) && ~all(coeff(:) == 0)) || isa(coeff, 'chebfun2') || isa(coeff, 'function_handle') )
        % There's nothing to do here.
        continue
    end
    
    % Define operators
    CC{k} = cell(length(coeff), 2);
    Sx = util.convertmat(n, dx, 1);
    Sy = util.convertmat(n, dy, 1);
    Dx = (2/diff(dom(1:2)))^dx * util.diffmat(n, dx);
    Dy = (2/diff(dom(3:4)))^dy * util.diffmat(n, dy);
    
    if ( isnumeric(coeff) )
        % Constant coefficient
        CC{k}{1,1} = coeff .* Sy * Dy;
        CC{k}{1,2} =          Sx * Dx;
    else
        % Variable coefficient
        % Compute a CDR decomposition:
        [C, D, R] = cdr(coeff, n, dom);
        % Make a multiplication operator for each slice
        rk = size(D, 1);
        CC{k} = cell(rk, 2);
        for r = 1:rk
            Mx = util.multmat( n, sign(D(r,r)) * sqrt(abs(D(r,r))) * R(:,r), dx );
            My = util.multmat( n,                sqrt(abs(D(r,r))) * C(:,r), dy );
            CC{k}{r,1} = Sy * My * Dy;
            CC{k}{r,2} = Sx * Mx * Dx;
        end
    end
end

% Flatten into a rank x 2 cell array
CC = cat(1, CC{:});

end

function [C, D, R] = cdr(f, n, dom)
%CDR   Compute a low rank decomposition of a bivariate function.
%   [C, D, R] = cdr(F, N, DOM) produces a K x K diagonal matrix D and N x K
%   matrices C and R of size such that F(x,y) = C(y,:) * D * R(x,:)'. Here,
%   K is the computed rank of the function F over the domain DOM, and each
%   column of C and R contains the Chebyshev coefficients of a univariate
%   function.

if ( isa(f, 'chebfun2') )
    [C, D, R] = cdr(f);
    C = chebcoeffs(C, n);
    R = chebcoeffs(R, n);
else
    [xx, yy] = util.chebpts2(n, n, dom);
    A = f(xx,yy);

    % Old way based on SVD: (This gives bad roundoff errors.)
    %   [C, D, R] = svd(A);
    %   r = rank(D);

    % GE with complete pivoting: (This is more accurate.)
    scl = max( abs( A(:) ) );
    C = []; D = []; R = [];
    for r = 1:n
        [mx, idx] = max( abs( A(:) ) ); % Complete pivoting
        if ( mx/scl < 100*eps )
            r = r - 1;
            break
        end
        [j, k] = ind2sub( size(A), idx );
        C = [C A(:,k)];
        D(r,r) = 1/A(j,k);
        R = [R A(j,:).'];
        A = A - A(:,k)*A(j,:)/A(j,k);
    end

    D = D(1:r,1:r);
    C = util.vals2coeffs(C(:,1:r));
    R = util.vals2coeffs(R(:,1:r));

    C(abs(C*sqrt(D)) < eps * max(abs(C*sqrt(D)))) = 0;
    R(abs(R*sqrt(D)) < eps * max(abs(R*sqrt(D)))) = 0;
end

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
%ZERODOF  Eliminate some degrees of freedom in the matrix equation can be removed.

if ( nargin < 6)
    trans = false;
end

tol = 100*eps;
n = size(C1,1);
perm = [2 1 3];
if ( trans )
    perm = [1 2 3];
end

for ii = 1:size(B, 1) % For each boundary condition, zero a column.
    C1ii = C1(:,ii); % Constant required to zero entry out.
    if ( ~any( abs(C1ii) > tol ) ), continue, end
    C1 = C1 - C1ii*sparse(B(ii,:));
    Gii = permute(G(ii,:,:), [2 3 1]);
    C2Gii = C2*Gii;
    R = repelem(full(C1ii), n, 1) .* repmat(C2Gii, n, 1);
    R = reshape(R, n, n, 4*n);
    R = permute(R, perm);
    E = E - R;
end

C1(abs(C1) < tol) = 0;

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

function P = compatibleProjection(p)
%COMPATIBLEPROJECTION Construct an orthogonal projection onto compatible
% boundary conditions.
%
%   P is a 4p x 4p matrix that imposes continuity at the corners.

% Order is: [lbc; rbc; dbc; ubc]

% Top left corner
A = zeros(4,4*p);
A(1,1:p) = ones(1,p);
A(1,3*p+1:end) = -(-1).^(0:p-1);

% Bottom left corner
A(2,1:p) = (-1).^(0:p-1);
A(2,2*p+1:3*p) = -(-1).^(0:p-1);

% Top right corner
A(3,p+1:2*p) = ones(1,p);
A(3,3*p+1:end) = -ones(1,p);

% Bottom right corner
A(4,p+1:2*p) = (-1).^(0:p-1);
A(4,2*p+1:3*p) = -ones(1,p);

% Build the projection:
VV = null(A);
P = VV * VV';

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% UTILITY CODES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

function [S22, Ainv] = mysolve(A, BC, rhs, p, pref)
% This involves solving a O(p^2) x O(p^2) almost-banded-block-banded system
% O(p) times, which we can do in O(p^4) using Schur complements / Woodbury.

if ( nargin < 5 )
    pref = ultraSEM.Pref();
end

switch pref.solver
    case '\'
        % Simply use backslash:
        S22 = A \ [BC, rhs];
        Ainv = @(u) A\u;
    case 'woodbury'
        % Woodbury formula:
        % The number of dense rows scales with the bandwidth of A
        m = bandwidth(A, 'lower');
        S22 = schurSolve(A, [BC, rhs], m);
        Ainv = @(u) schurSolve(A, u, m);
    case 'LU'
        % Do sparse LU by hand so we can store L U factors:
        P = symrcm(A);
        [~, Pinv] = sort(P);
        [L, U, p] = lu(A(P,P), 'vector');
        rowPermute = @(v,p) v(p,:);
        Ainv = @(b) rowPermute((U\(L\b(P(p),:))), Pinv);
        S22 = Ainv([BC, rhs]);
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

AA = A(i2,i2);
bb = [b(i2,:), full(A(i2,i1))];
if ( nnz(AA)/numel(AA) > .7 )
    % A is dense. Bail and do full direct solve.
    x = full(A)\b;
    return
end

if ( doRowScaling )
    % Row scaling to improve accuracy
    s = 1./ max(1, max(abs(AA), [], 2) );
    AA = bsxfun(@times, s, AA);
    bb = s.*bb;
end

% Force banded solver:
parms = spparms;
spparms('bandden', 0);
c = AA\bb;
spparms(parms);

x = (full(A(i1,i1)) - full(A(i1,i2))*c(:,i3)) \ (b(i1,:) - full(A(i1,i2))*c(:,1:nb));
y = c(:,1:nb) - c(:,i3)*x;
x = [x ; y];

end

function uniformOp = isUniformOp(op, dom)
% Operator is 'uniform' if it is constant on each patch, and each patch has
% the same constant Jacobian (i.e., rectangular patches of uniform size).

uniformOp = false;
if ( isa(dom, 'ultraSEM.Rect') )
    % TODO: This breaks encapsulation. Ignore it for now.
    dom = rectVertices(dom);
end
if ( isnumeric(dom) )
    [~, isConstant] = feval(op, dom(1,1), dom(1,2));
    if ( isConstant && ~any(diff(diff(dom(:,1:2),1,2))) && ...
            ~any(diff(diff(dom(:,3:4),1,2))) )
        uniformOp = true;
    end
end
end