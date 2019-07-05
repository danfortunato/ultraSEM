classdef ultraSEM < handle
%ULTRASEM class.
%   ULTRASEM(DOM, OP, RHS) constructs an ULTRASEM object from the
%   ultraSEMDomain DOM and the partial differential operator (PDO) OP. Here
%   OP may be either an ultraSEMPDO object or a cell array containing the
%   information which would otherwise be passed to the ultraSEMPDO
%   constructor. (See the documentation of the ultraSEMPDO class for
%   details.)
%
%   ULTRASEM(DOM, OP) assumes the problem is homogeneous (i.e., RHS = 0).
%
%   ULTRASEM(DOM, OP, RHS, N) uses an N x N discretization to solve on each
%   patch of the domain DOM.
%
%   The full sequence for solving a problem using an ultraSEM object S is
%       S = ultraSEM(DOM, OP, RHS)
%       build(S)
%       sol = S\bc % or sol = solve(S, bc)
%   However, this can be abbreviated to
%       S = ultraSEM(DOM, OP, RHS)
%       sol = S\bc % or sol = solve(S, bc)
%   The advantage of the former is that the initialization and build phases
%   (in which most of the computation work is performed) can be reused for
%   different boundary conditions.
%
% Example:
%   S = ultraSEM.alphabet('S') + 2i;     % Form an 'S' domain.
%   U = ultraSEM.alphabet('U') + 2;      % Form a 'U' domain.
%   SU = S & U;                          % Combine the two.
%   SU = refine(SU, 1);                  % Refine the grid.
%   op = ultraSEMPDO(1, 0, @(x,y) 50*y); % Construct a PDO: lap(u) + 50*y*u.
%   rhs = -1;                            % RHS (scalar or function handle).
%   bc = 0;                              % Boundary condition (as above).
%   L = ultraSEM(SU, op, rhs);           % Construct the ultraSEM object.
%   sol = L\bc;                          % Solve.
%   plot(sol)                            % Plot.
%
% The typical user workflow may presented as:
%
%  [Specify domain   ]                 [Specify PDO   ]
%  [as ultraSEMDomain]                 [as ultraSEMPDO]
%          \                             /
%           \                           /
%             -> [Construct ultraSEM]<-
%                         |
%                      (build)       <-- optional
%                         |
%                       solve
%                         |
%                         v
%                   [ultraSEMSol]

% Copyright 2018 by Nick Hale and Dan Fortunato.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DEVELOPER NOTES:
% * Note that ULTRASEM is a handle class, so modifications to its
%   properties made in methods are altered globally.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties ( Access = public )

        domain  = [] % Domain of the ultraSEM object. (An ultraSEMDomain)
        patches = {} % Cell array of ultraSEMPatches.

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods

        function obj = ultraSEM(dom, op, rhs, n)
        %ULTRASEM  Constructor for the ultraSEM class.
        %   Documentation at the head of ultraSEM.m file.

            if ( nargin == 0 )
                % Construct empty object:
                return
            end

            % Assign the domain:
            obj.domain = dom;

            if ( nargin < 3 )
                % Default to homogeneous problem:
                rhs = 0;
            end
            if ( nargin < 4 )
                % Set discretization size:
                n = [];
            end

            % Initialize patches:
            obj.initialize(op, rhs, n);

        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% METHODS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods ( Access = public, Static = false )
        
        function coeffs = bc2coeffs(S, bc)
        %BC2COEFFS   Convert boundary conditions to coeffs.

            if ( numel(S.patches) ~= 1 )
                error('ULTRASEM:ULTRASEM:bc2coeffs:notBuilt', ...
                    'The ultraSEM object %s must be built before creating BCs.', ...
                    inputname(1));
            end

            % Get boundary points:
            xy = S.patches{1}.xy;

            % Initialize boundary data. We should have n coefficients for
            % each boundary.
            coeffs = cell(size(xy));
            if ( ~isnumeric(bc) )
                % Evaluate the BC if given a function handle:
                for k = 1:size(coeffs, 1)
                    vals = feval(bc, xy{k}(:,1), xy{k}(:,2));
                    % Convert from values to coeffs:
                    coeffs{k} = util.vals2coeffs(vals);
                end
            elseif ( isscalar(bc) )
                % Convert a scalar to coeffs:
                for k = 1:size(coeffs, 1)
                    n = size(xy{k}, 1);
                    coeffs{k} = [bc ; zeros(n-1, 1)];
                end
            else
                error('ULTRASEM:ULTRASEM:bc2coeffs:badBC', ...
                    'Cannot evaluate boundary data.');
            end

            % CAT() is 10x faster than CELL2MAT().
            coeffs = cat(1, coeffs{:});

        end

        function build(S)
        %BUILD   Build the merge tree of patches in an ultraSEM object.
        %   BUILD(S) will construct the global solution operator from the
        %   local operators in S.patches as in the order defined in the
        %   merge tree S.domain.mergeIdx. After building, S will have a
        %   single patch.
        %
        %   If S has not yet been initialized, then an error is thrown.
        %
        % See also INITIALIZE, SOLVE.
        
            if ( isa(S.patches, 'ultraSEM') )
                % Recurse down and build lower level patches
                for k = 1:numel(S.patches)
                    build(S.patches(k));
                end
                % Concatenate (since S already contains domain info.)
                S.patches = vertcat(S.patches.patches);
            end

            if ( isempty(S.patches) || isempty(S.patches{1}.S) )
                error('ULTRASEM:ULTRASEM:build:notInitialized', ...
                    '%f has not yet been initialized.', inputname(1))
            end

            % Build the patches:
            S.patches = build(S.domain, S.patches);

        end

        function initialize(S, varargin)
        %INITIALIZE   Initialize an ultraSEM with a PDE and solve local subproblems.
        %   INITIALIZE(S, OP, RHS) will initialize the ultraSEM object S
        %   with the PDE defined by OP and the righthand side RHS. The PDE
        %   will be solved on each of the subpatches of S.domain.
        %
        %   INITIALIZE(S, OP, RHS, N) uses an N x N discretization on each
        %   patch.
        %
        %   The full sequence for solving a problem using an ultraSEM
        %   object S is
        %       initialize(S, OP, RHS)
        %       build(S)
        %       sol = S\bc % or sol = solve(S, bc)
        %
        % See also BUILD, SOLVE.

            % Initialize all leaf patches:
            
            D = S.domain;
            numD = size(D.domain, 1);
            if ( isa(D.domain, 'ultraSEMDomain') && numD > 1)
                S_sub(numD,1) = ultraSEM;
                for k = 1:numD
                    S_sub(k).domain = D.domain(k);
                    initialize(S_sub(k), varargin{:});
                end
                S.patches = S_sub;
            else
                S.patches = ultraSEMLeaf.initialize(D, varargin{:});
            end
            
        end

        function h = merge(f, g)
        %MERGE   Merge two ultraSEM objects.
        %   H = MERGE(F, G) will merge the two ultraSEM objects F and G. If
        %   F and G have already been initialized and built, then H will be
        %   also.
        %
        % See also INITIALIZE, BUILD.

            % Merge the domains:
            h = ultraSEM();
            h.domain = merge(f.domain, g.domain);

            % Merge the patches:
            if ( numel(f.patches) == 1 && numel(g.patches) == 1)
                % f and g have already been built.
                h.patches{1} = merge(f.patches{1}, g.patches{1});
            else
                h.patches = [f.patches ; g.patches];
                h.initialize;
            end

        end

        function varargout = mldivide(varargin)
        %MLDIVIDE   Solve an ultraSEM system.
        %   MLDIVIDE is a shorthand for SOLVE().
        %
        % See also SOLVE.

            % MLDIVIDE() is simply a wrapper for SOLVE():
            [varargout{1:nargout}] = solve(varargin{:});

        end

        function sol = solve(S, bc)
        %SOLVE   Solve an ultraSEM system.
        %   SOL = SOLVE(S, BC) returns an ultraSEMSol object SOL
        %   representing the solution to the PDE, subject to the boundary
        %   conditions specified by the function handle BC, contained in
        %   the ultraSEM object S. If S has not yet been built (see
        %   BUILD()) then SOLVE() will build it. If S has not yet been
        %   initialized (see INITIALIZE()) then an error is thrown.
        %
        %   The full sequence for solving a problem using an ultraSEM
        %   object S is
        %       initialize(S, OP, RHS)
        %       build(S)
        %       sol = S\bc % or sol = solve(S, bc)
        %
        % See also BUILD, INITIALIZE.

            % Build the ultraSEM if required:
            if ( numel(S.patches) > 1 ), build(S); end

            if ( numel(S.patches) == 0 )
                error('ULTRASEM:ULTRASEM:solve:notInitialized', ...
                    'The ultraSEM object `%s` has not been initialized.', ...
                    inputname(1));
            end

            % Build the boundary conditions:
            coeffs = bc2coeffs(S, bc);

            % Solve the patch object:
            sol = solve(S.patches{1}, coeffs);

        end

        function S = subsasgn(S, index, val)
        %SUBSASGN   Subscripted assignment to an ultraSEM object.
        %   S.PROP = VAL allows property assignment for S.DOMAIN and
        %   S.PATCHES.
        %
        %   S.RHS = F is equivalent to UPDATERHS(S, F).
        %
        %   S() and S{} are not supported.
        %
        % See also UPDATERHS.

            switch index(1).type
                case '.'

                    if ( strcmpi(index(1).subs, 'rhs') )
                        S = updateRHS(S, val);
                    else
                        S = builtin('subsasgn', S, index, val);
                    end

                otherwise
                    error('ULTRASEM:ULTRASEM:subsasgn:UnexpectedType', ...
                        ['??? Unexpected index.type of ' index(1).type]);
            end

        end

        function S = updateRHS(S, f)
        %UPDATERHS   Update the RHS function in a ultraSEM object.
        %   UPDATERHS(S, F) or S.RHS = F updates the RHS function F in the
        %   ultraSEM object S. This is useful as only a small fraction of
        %   the initialization phase needs to be repeated.
        %
        % Example:
        %   S = ultraSEM.alphabet('S');
        %   T = refine(S, 2);
        %   op = {1, 0, @(x,y) y};
        %   tic
        %       S = ultraSEM(T, op, -1);
        %       S.build;
        %       sol1 = S\0;
        %   t1 = toc;               % t1 = 0.866088
        %   tic
        %       S.rhs = @(x,y) sin(pi*x.*y);
        %       sol2 = S\0;
        %   t2 = toc;               % t2 = 0.277288
        %
        % See also BUILD.

            for k = 1:numel(S.patches)
                S.patches{k} = updateRHS(S.patches{k}, f);
            end

        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STATIC METHODS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods ( Access = public, Static = true )

        function varargout = alphabet(varargin)
        %ALPHABET   Construct an ultraSEMDomain alphabet domain.
        %   ULTRASEM.ALPHABET(STR), where STR is a single string character,
        %   returns an ultraSEMDomain object in the shape of STR.
            [varargout{1:nargout}] = ultraSEMDomain.alphabet(varargin{:});
        end

        function T = delaunay(v, P)
        %DELAUNAY   Delaunay triangularization.
        %   T = ULTRASEM.DELAUNAY(X), where is an N x 2 matrix, returns an
        %   ultraSEMDomain representation of the Delaunay triangulation
        %   equivalent to DELAUNAY(X), but with each triangle subdivided
        %   into three further kites with their intersection at the
        %   barycenter.
        %
        %   T = ULTRASEM.DELAUNAY(X_BDY, X_INT), where X_INT is an N x 2
        %   matrix and X_BDY is an M x 2 matrix, creates a constrained
        %   Delaunay triangulation, equivalent to
        %   delaunayTriangulation(X_BDY, X_INT).
        %
        % See also TRIANGLE, QUAD.

            % Compute triangulariztation:
            if ( nargin < 2 )
                dt = delaunayTriangulation(v);
                list = dt.ConnectivityList;
            else
                dt = delaunayTriangulation(v, P);
                IO = isInterior(dt);
                list = dt(IO,:);
            end
            pts = dt.Points(:,:);

            %% Build domain:
            nt = size(list,1);
            T = [];
            for k = 1:nt
                v = pts(list(k,:),:);
                T = T & ultraSEM.triangle(v);
            end

        end

        function varargout = rectangle(varargin)
        %RECTANGLE   Construct a rectangular ultraSEMDomain.
        %   ULTRASEM.RECTANGLE([A B C D]) returns a rectangluar
        %   ultraSEMDomain with vertices (A,C), (A,D), (B,D), (B,C) with a
        %   single patch.
        %
        %   ULTRASEM.RECTANGLE([A B C D], M, N) returns the same rectangle,
        %   but subdivided into an MxN grid of equally sized patches.
        %   Similarly, ULTRASEM.RECTANGLE([A B C D], M) assumes N = M.
            [varargout{1:nargout}] = ultraSEMDomain.rectangle(varargin{:});
        end

        function varargout = ultraSEMQuad(varargin)
        %QUAD Construct an ultraSEMDomain quadrilateral domain.
        %   ultraSEM.ultraSEMQuad(V) constructs a quadrilateral with vertices V.
            [varargout{1:nargout}] = ultraSEMDomain.ultraSEMQuad(varargin{:});
        end

        function varargout = triangle(varargin)
        %TRIANGLE  Construct an ultraSEMDomain triangle domain.
        %   ultraSEM.triangle(V) constructs a triangle with vertices V.
            [varargout{1:nargout}] = ultraSEMDomain.triangle(varargin{:});
        end

        function varargout = polygon(varargin)
        %POLYGON  Construct an ultraSEMDomain convex polygon domain.
        %   ultraSEM.polygon(V) constructs a convex polygon with vertices
        %   V.
            [varargout{1:nargout}] = ultraSEMDomain.polygon(varargin{:});
        end

        function varargout = trimesh(varargin)
        %TRIMESH  Construct an ultraSEMDomain triangulated mesh.
        %   ultraSEM.trimesh(P, T) constructs a triangulated mesh.
            [varargout{1:nargout}] = ultraSEMDomain.trimesh(varargin{:});
        end

        % Run the test suite.
        pass = test();
        
        % Run the test suite.
        duration = bench();

        % Run the GUI.
        varargout = gui();

    end

end
