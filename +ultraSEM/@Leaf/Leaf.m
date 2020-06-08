classdef Leaf < ultraSEM.Patch
%ULTRASEM.LEAF   Leaf subclass of a patch (where subproblems are solved).
%   P = ULTRASEM.LEAF(DOMAIN, S, D2N, EDGES, AINV) creates an ULTRASEM.LEAF
%   object P and assigns each of the inputs to their associated properties
%   in P.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DEVELOPER NOTES:
% * Note: We assume that an ULTRASEM.LEAF has four sides and that its
%   boundary nodes XY are stored in the order "left", "right", "down",
%   "up".
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties ( Access = public )

        Ainv     % Local (homogeneous BC) solution operator.
        normal_d % Normal derivative operator.
                 % (These are stored so the RHS can be efficiently updated)

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods

        function P = Leaf(dom, S, D2N, edges, Ainv, normal_d)

            % Construct empty patch:
            if ( nargin == 0 )
                return
            end

            % Assign domain and operators:
            P.domain = dom;            % Domain.
            P.S = S;                   % Solution operator.
            P.D2N = D2N;               % Dirichlet-to-Neumann map.
            P.edges = edges;           % Edges.
            if ( nargin > 4 )
                P.Ainv = Ainv;         % Local solution operator.
            end
            if ( nargin > 5 )
                P.normal_d = normal_d; % Normal derivative operator.
            end

            % The D2N map needs to be scaled on each side (e.g. when being
            % merged) to account for the Jacobian scaling which has been
            % factored out of the coordinate derivative maps. This scaling
            % is not known until the merge stage, as it depends on the
            % scaling of the neighboring patch.
            if ( isa(dom, 'ultraSEM.Rect') )
                % Normal derivatives are computed differently for Rect.
                % Don't scale.
                P.D2N_scl = {1, 1, 1, 1}.';
            else
                P.D2N_scl = cell(4,1);
                P.D2N_scl{1} = @(y) dom.det(-1, y); % Left
                P.D2N_scl{2} = @(y) dom.det(1, y);  % Right
                P.D2N_scl{3} = @(x) dom.det(x, -1); % Down
                P.D2N_scl{4} = @(x) dom.det(x, 1);  % Up
            end

        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STATIC METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods ( Static )

        % Initialize an array of ULTRASEM.LEAF objects.
        P = initialize(dom, op, varargin);

    end

end
