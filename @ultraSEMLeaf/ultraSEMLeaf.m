classdef ultraSEMLeaf < ultraSEMPatch
%ULTRASEMLEAF  Leaf subclass of an ultraSEMPatch (where subproblems are solved).

% Copyright 2018 by Nick Hale and Dan Fortunato.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DEVELOPER NOTES:
% * Note: We assume that an ULTRASEMLEAF has four sides and that its
%   boundary nodes XY are stored in the order "left", "right", "down",
%   "up".
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %#ok<*PROP>

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties ( Access = public )

        n  % Discretization size.
        op % Description of PDO (stored so the RHS can be efficiently updated).

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods

        function P = ultraSEMLeaf(domain, S, D2N, xy, n, op)
        %ULTRASEMLEAF   Class constructor for the @ultraSEMLeaf class.
        %   P = ultraSEMLeaf(DOMAIN, S, D2N, XY, N, OP) assigns each of
        %   the inputs to their associated properties in the ultraSEMLeaf
        %   object OBJ.

            % Construct empty patch:
            if ( nargin == 0 )
                return
            end

            % Assign domain and operators:
            P.domain = domain;    % Domain.
            P.S = S;              % Solution operator.
            P.D2N = D2N;          % Dirichlet-to-Neumann map.
            P.xy = xy;            % Boundary nodes.
            P.n = n;              % Discretization size in x.
            P.op = op;            % PDO (in cell form).

        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% METHODS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods ( Access = public, Static = false )

        function [u, d] = solve(P, bc)
        %SOLVE   Solve a leaf patch.
        %   SOL = SOLVE(P, BC) returns an ultraSEMSol object representing
        %   the PDE solution on the patch P with Dirichlet boundary data
        %   BC.
        %
        %   [U, D] = SOLVE(P, BC) returns instead a cell array containing
        %   the solution coefficients U and a vector D containing the
        %   domain of the patch.

            % Extract the domain from the patch:
            d = P.domain;

            % Evaluate the solution operator for the patch:
            u = P.S * [bc ; 1]; % The 1 accounts for the particular part.
            u = reshape(u, P.n, P.n);

            % Return cell output for consistency with
            % ultraSEMParent/solve():
            u = {u};

            if ( nargout == 1 )
                % If a single output is requested, return an ultraSEMSol
                % object.
                u = ultraSEMSol(u, d);
            end

        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STATIC METHODS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods ( Access = public, Static = true )

        % Initialize an array of ultraSEMLeaf objects.
        P = initialize(dom, op, rhs, n);

        function P = updateRHS(P, rhs)
        %UPDATERHS   Update RHS of a spceddLeaf object.
        %   P = UPDATERHS(P, F) replaces the existing RHS of an initialised spceddLeaf
        %   object P with that given in F, which may be a constant or a function handle.
        %   Since only one subproblem needs to be solved on each patch in the update
        %   process (rather than O(n) in the original intitialisation) this can lead to
        %   considerable performance gains when solving for multiple RHSs.
        %
        % See also SPECDDLEAF.INITIALISE.

            error('ULTRASEM:ULTRASEMLEAF:updateRHS:notImplemented', ...
                'Updating the RHS of a leaf has not been implemented yet.');

        end

    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% METHODS IMPLEMENTED IN THIS FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
