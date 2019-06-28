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

        Ainv  % Local (homogeneous BC) solution operator operator.
              % (Stored so the RHS can be efficiently updated)
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods

        function P = ultraSEMLeaf(dom, S, D2N, xy, Ainv)
        %ULTRASEMLEAF   Class constructor for the @ultraSEMLeaf class.
        %   P = ultraSEMLeaf(DOMAIN, S, D2N, XY) assigns each of the inputs to
        %   their associated properties in the ultraSEMLeaf object OBJ.

            % Construct empty patch:
            if ( nargin == 0 )
                return
            end

            % Assign domain and operators:
            P.domain = dom;       % Domain.
            P.S = S;              % Solution operator.
            P.D2N = D2N;          % Dirichlet-to-Neumann map.
            P.xy = xy;            % Boundary nodes.
            if ( nargin > 4 )
                P.Ainv = Ainv;    % Local solution operator.
            end

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
            u = reshape(u, P.n*[1,1]);

            % Return cell output for consistency with
            % ultraSEMParent/solve():
            u = {u};

            if ( nargout == 1 )
                % If a single output is requested, return an ultraSEMSol
                % object.
                u = ultraSEMSol(u, d);
            end

        end
        
        function out = n(P)
            out = size(P.D2N,1)/4;
        end

    end
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STATIC METHODS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods ( Access = public, Static = true )

        % Initialize an array of ultraSEMLeaf objects.
        P = initialize(dom, op, rhs, n);

    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% METHODS IMPLEMENTED IN THIS FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
