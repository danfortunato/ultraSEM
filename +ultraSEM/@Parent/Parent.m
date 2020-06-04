classdef Parent < ultraSEM.Patch
%ULTRASEM.PARENT   Parent subclass of a patch.
%   P = ULTRASEM.PARENT(DOMAIN, S, D2N, EDGES, CHILD1, CHILD2, IDX1, IDX2)
%   creates an ULTRASEM.PARENT object P and assigns each of the inputs to
%   their associated properties in P.

% Copyright 2018 by Nick Hale and Dan Fortunato.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties ( Access = public )

        child1 = [] % Child patch
        child2 = [] % Child patch
        idx1        % How p.xy relates to p.child1.xy
        idx2        % How p.xy relates to p.child2.xy
        l2g1        % Local-to-global map from child1 to interface
        l2g2        % Local-to-global map from child2 to interface
        scl1
        scl2
        dA          % Decomposition of interface linear system

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods

        function P = Parent(domain, S, D2N, dA, edges, scl, child1, child2, idx1, idx2, l2g1, l2g2, scl1, scl2)

            % Construct empty patch:
            if ( nargin == 0 )
                return
            end

            % Assign domain and operators:
            P.domain = domain;
            P.S = S;
            P.D2N = D2N;
            P.dA = dA;
            P.edges = edges;
            P.D2N_scl = scl;

            % Assign children:
            P.child1 = child1;
            P.child2 = child2;
            P.idx1 = idx1;
            P.idx2 = idx2;
            P.l2g1 = l2g1;
            P.l2g2 = l2g2;
            P.scl1 = scl1;
            P.scl2 = scl2;

        end

    end

end