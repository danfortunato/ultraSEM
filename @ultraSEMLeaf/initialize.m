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
corners = [n-1 n];
oneBdyDOF = n-2;
numBdyDOF = 4*oneBdyDOF;
numIntDOF = n^2-4;
numRHSDOF = numIntDOF - numBdyDOF;

ii = 1:n^2; ii([n^2-n-1 n^2-n n^2-1 n^2]) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%% DEFINE BOUNDARY NODES %%%%%%%%%%%%%%%%%%%%%%%%%%%%
[XX, YY] = chebpts2(n);          % Chebyshev points and grid.
leftIdx  = sub2ind([n n], (1:n).', ones(n,1));
rightIdx = sub2ind([n n], (1:n).', n*ones(n,1));
downIdx  = sub2ind([n n], ones(n,1), (1:n).');
upIdx    = sub2ind([n n], n*ones(n,1), (1:n).');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%% DEFINE BOUNDARY CONDITIONS %%%%%%%%%%%%%%%%%%%%%%%%%
I = speye(n);

% Construct Dirichlet conditions along the four edges:
lbc = kron( (-1).^(0:n-1), I );
rbc = kron( ones(1,n), I );
dbc = kron( I, (-1).^(0:n-1) );
ubc = kron( I, ones(1,n) );
% Remove the corner conditions:
lbc(corners,:) = [];
rbc(corners,:) = [];
dbc(corners,:) = [];
ubc(corners,:) = [];
bcrows = [ lbc ; rbc ; dbc ; ubc ];

% Construct normal derivatives conditions along the four edges:
lbc_d = kron( (-1).^(0:n-1).*(0:n-1).^2, I );
rbc_d = kron( ones(1,n).*(0:n-1).^2, I );
dbc_d = kron( I, (-1).^(0:n-1).*(0:n-1).^2 );
ubc_d = kron( I, ones(1,n).*(0:n-1).^2 );
% Again, we remove the corner conditions. Normal derivatives are not
% well-defined at corners anyway:
lbc_d(corners,:) = [];
rbc_d(corners,:) = [];
dbc_d(corners,:) = [];
ubc_d(corners,:) = [];
bcrows_d = [ lbc_d ; rbc_d ; dbc_d ; ubc_d ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%% DEFINE OPERATORS %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Scaling (for all patches):
domk = dom(1,:);
domx = domk(1:2);
domy = domk(3:4);
sclx = 2/diff(domx);
scly = 2/diff(domy);

% Define operators (on [-1,1]^2):
D   = ultraS.diffmat( n, 1 );
D2  = ultraS.diffmat( n, 2 );       D2  = D2(1:end-2,:);
S01 = ultraS.convertmat( n, 0, 1 ); S01 = S01(1:end-2,:);
S1  = ultraS.convertmat( n, 1, 1 );
S1D = S1*D; S1D = S1D(1:end-2,:);
Lxx = kron(D2,  S01);
Lyy = kron(S01, D2);
Lxy = kron(S1D, S1D);
Lx  = kron(S1D, S01);
Ly  = kron(S01, S1D);
L0  = kron(S01, S01);

% Construct the differential operator:
A = (sclx^2*op.dxx)*Lxx + (scly^2*op.dyy)*Lyy + (sclx.*scly*op.dxy)*Lxy + ...
       (sclx*op.dx)*Lx  +    (scly*op.dy)*Ly  + op.b*L0;
A = [ bcrows(:,ii) ; A(:,ii) ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% CONSTANT RHS? %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define scalar RHSs:
if ( isnumeric(rhs) )
    constantRHS = true;
    if ( isscalar(rhs) )
        rhs_eval = [zeros(numBdyDOF,1); rhs; zeros(numRHSDOF-1,1)];
    else
        % The user gave us coeffs:
        rhs_eval = [zeros(numBdyDOF,1); rhs(1:numRHSDOF)];
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
BC = eye(numIntDOF, numBdyDOF); % Encode every possible BC.
Si = A \ [BC, rhs_eval(:,1)];
% Append boundary points to solution operator:
S = zeros(n^2, size(Si,2));     % Set the missing 4 highest modes to zero
S(ii,:) = Si;                   % (due to the corners).

% Dirichlet-to-Neumann map:
D2N = bcrows_d * S;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize
P = cell(numPatches, 1);

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
        % Drop the four highest modes (due to the corners):
        rhs_eval(:,k) = [zeros(numBdyDOF,1); coeffs(1:numRHSDOF).'];
    end

    % Assemble the patch:
    P{k} = ultraSEMLeaf(dom(k,:), S, D2N, xy, n, op);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Append particular parts:
if ( ~constantRHS )
    Si = A \ rhs_eval;
    S = zeros(n^2, size(Si,2));
    S(ii,:) = Si;
    D2N = bcrows_d * S;
    for k = 1:numPatches
        P{k}.S(:,end) = S(:,k);
        P{k}.D2N(:,end) = P{k}.D2N(:,end) + D2N(:,k);
    end
end

end
