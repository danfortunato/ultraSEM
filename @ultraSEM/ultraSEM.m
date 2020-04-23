classdef ultraSEM < handle
%ULTRASEM class.
%   ULTRASEM(DOM, OP, RHS) constructs an ULTRASEM object from the
%   ULTRASEM.DOMAIN DOM and the partial differential operator (PDO) OP.
%   Here OP may be either an ULTRASEM.PDO object or a cell array containing
%   the information which would otherwise be passed to the ULTRASEM.PDO
%   constructor. (See the documentation of the ULTRASEM.PDO class for
%   details.)
%
%   ULTRASEM(DOM, OP) assumes the problem is homogeneous (i.e., RHS = 0).
%
%   ULTRASEM(DOM, OP, RHS, N) uses an N x N discretization to solve on each
%   patch of the domain DOM.
%
%   ULTRASEM(..., PREF) uses the preferences specified in the ULTRASEM.PREF
%   object PREF. (See ULTRASEM.PREF for details on the various preference
%   options and their defaults.)
%
%   The full sequence for solving a problem using an ULTRASEM object S is:
%
%       S = ultraSEM(DOM, OP, RHS)
%       build(S)
%       sol = S\bc % or sol = solve(S, bc)
%
%   However, this can be abbreviated to:
%
%       S = ultraSEM(DOM, OP, RHS)
%       sol = S\bc % or sol = solve(S, bc)
%
%   The advantage of the former is that the initialization and build phases
%   (in which most of the computation work is performed) can be reused for
%   different boundary conditions.
%
% Example:
%
%   S = ultraSEM.alphabet('S') + 2i;      % Form an 'S' domain.
%   U = ultraSEM.alphabet('U') + 2;       % Form a 'U' domain.
%   SU = S & U;                           % Combine the two.
%   SU = refine(SU, 1);                   % Refine the grid.
%   op = ultraSEM.PDO(1, 0, @(x,y) 50*y); % Construct a PDO: lap(u) + 50*y*u.
%   rhs = -1;                             % RHS (scalar or function handle).
%   bc = 0;                               % Boundary condition (as above).
%   L = ultraSEM(SU, op, rhs);            % Construct the ULTRASEM object.
%   sol = L\bc;                           % Solve.
%   plot(sol)                             % Plot.
%
% The typical user workflow may presented as:
%
%  [Specify domain    ]               [Specify PDO    ]
%  [as ULTRASEM.DOMAIN]               [as ULTRASEM.PDO]
%          \                             /
%           \                           /
%             -> [Construct ULTRASEM]<-
%                         |
%                      (build)       <-- optional
%                         |
%                       solve
%                         |
%                         v
%                   [ULTRASEM.SOL]

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

        domain  = [] % Domain of the ultraSEM object. (An ultraSEM.Domain)
        patches = {} % Cell array of ultraSEM.Patches.

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods

        function obj = ultraSEM(dom, op, varargin)
        %ULTRASEM   Constructor for the ULTRASEM class.
        %   Documentation at the head of ultraSEM.m file.

            if ( nargin == 0 )
                % Construct empty object:
                return
            end

            % Assign the domain:
            obj.domain = dom; % TODO: Validate domain.

            % Set defaults for unspecified arguments:
            rhs = 0; % Default to homogeneous problem.
            n = [];
            pref = ultraSEM.Pref();

            if ( nargin == 3 )
                if ( isa(varargin{1}, 'ultraSEM.Pref') )
                    pref = varargin{1}; % ULTRASEM(DOM, OP, PREF)
                else
                    rhs = varargin{1};  % ULTRASEM(DOM, OP, RHS)
                end
            elseif ( nargin == 4 )
                rhs = varargin{1};
                if ( isa(varargin{2}, 'ultraSEM.Pref') )
                    pref = varargin{2}; % ULTRASEM(DOM, OP, RHS, PREF)
                else
                    n = varargin{2};    % ULTRASEM(DOM, OP, RHS, N)
                end
            elseif ( nargin == 5 )
                rhs = varargin{1};      % ULTRASEM(DOM, OP, RHS, N, PREF)
                n = varargin{2};
                pref = varargin{3};
            end

            % Initialize patches:
            obj.initialize(op, rhs, n, pref);

        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% PRIVATE METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods ( Access = private, Static = false )

        % Convert boundary conditions to Chebyshev coefficients.
        coeffs = bc2coeffs(S, bc);

        % Check to see if an ultraSEM has been initialized.
        out = isInitialized(S);
        
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STATIC METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods ( Access = public, Static = true )

        function varargout = alphabet(varargin)
        %ALPHABET   Construct an ULTRASEM.DOMAIN alphabet domain.
        %   ULTRASEM.ALPHABET(STR), where STR is a single string character,
        %   returns an ULTRASEM.DOMAIN object in the shape of STR.
            [varargout{1:nargout}] = ultraSEM.Domain.alphabet(varargin{:});
        end

        function varargout = rectangle(varargin)
        %RECTANGLE   Construct a rectangular ULTRASEM.DOMAIN.
        %   ULTRASEM.RECTANGLE([A B C D]) returns a rectangluar
        %   ULTRASEM.DOMAIN with vertices (A,C), (A,D), (B,D), (B,C) with a
        %   single patch.
        %
        %   ULTRASEM.RECTANGLE([A B C D], M, N) returns the same rectangle,
        %   but subdivided into an M x N grid of equally sized patches.
        %   Similarly, ULTRASEM.RECTANGLE([A B C D], M) assumes N = M.
            [varargout{1:nargout}] = ultraSEM.Domain.rectangle(varargin{:});
        end

        function varargout = quad(varargin)
        %QUAD Construct an ULTRASEM.DOMAIN quadrilateral domain.
        %   ULTRASEM.QUAD(V) constructs a quadrilateral with vertices V.
            [varargout{1:nargout}] = ultraSEM.Domain.quad(varargin{:});
        end

        function varargout = triangle(varargin)
        %TRIANGLE   Construct an ULTRASEM.DOMAIN triangle domain.
        %   ULTRASEM.TRIANGLE(V) constructs a triangle with vertices V by
        %   splitting a triangle into three quadrilaterals.
            pref = ultraSEM.Pref();
            if ( pref.splitTriangles )
                [varargout{1:nargout}] = ultraSEM.Domain.triangle(varargin{:});
            else
                [varargout{1:nargout}] = ultraSEM.Domain.duffy(varargin{:});
            end
        end

        function varargout = duffy(varargin)
        %DUFFY   Construct an ULTRASEM.DOMAIN triangle domain.
        %   ULTRASEM.DUFFY(V) constructs a triangle with vertices V based
        %   on the Duffy transformation.
            [varargout{1:nargout}] = ultraSEM.Domain.duffy(varargin{:});
        end

        function varargout = polygon(varargin)
        %POLYGON   Construct an ULTRASEM.DOMAIN convex polygon domain.
        %   ULTRASEM.POLYGON(V) constructs a convex polygon with vertices
        %   V.
            [varargout{1:nargout}] = ultraSEM.Domain.polygon(varargin{:});
        end

        function varargout = trimesh(varargin)
        %TRIMESH   Construct an ULTRASEM.DOMAIN triangulated mesh.
        %   ULTRASEM.TRIMESH(P, T) constructs a triangulated mesh.
            [varargout{1:nargout}] = ultraSEM.Domain.trimesh(varargin{:});
        end

        % Delaunay triangulation.
        T = delaunay(v, P);

        % Run the test suite.
        pass = test();
        
        % Run the benchmark suite.
        duration = bench();

        % Run the GUI.
        varargout = gui();

    end

end
