classdef ultraSEMParent < ultraSEMPatch

% Copyright 2018 by Nick Hale and Dan Fortunato.

    %#ok<*PROP>
    %#ok<*PROPLC>

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties ( Access = public )

        child1 = []     % Child patch
        child2 = []     % Child patch
        idx1            % How p.xy relates to p.child1.xy
        idx2            % How p.xy relates to p.child2.xy
        l2g1            % Local-to-global map from child1 to interface
        l2g2            % Local-to-global map from child2 to interface

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods

        function P = ultraSEMParent(domain, S, D2N, edges, child1, child2, idx1, idx2, l2g1, l2g2)
        %ULTRASEMPARENT   Class constructor for the @ultraSEMParent class.
        %   P = ultraSEMParent(DOMAIN, S, D2N, EDGES, CHILD1, CHILD2, IDX1, IDX2)
        %   assigns each of the inputs to their associated properties in
        %   the ultraSEMParent object P.

            % Construct empty patch:
            if ( nargin == 0 )
                return
            end

            % Assign domain and operators:
            P.domain = domain;
            P.S = S;
            P.D2N = D2N;
            P.edges = edges;

            % Assign children:
            P.child1 = child1;
            P.child2 = child2;
            P.idx1 = idx1;
            P.idx2 = idx2;
            P.l2g1 = l2g1;
            P.l2g2 = l2g2;

        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% METHODS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods ( Access = public, Static = false )

        function dof = dof(P)
        %DOF   Total number of degrees of freedom in an ultraSEMParent.
            dof = P.child1.dof + P.child2.dof;
        end

        function [u, d] = solve(P, bc)
        %SOLVE   Solve a patch.
        %   SOL = SOLVE(P, BC) returns an ultraSEMSol object representing
        %   the PDE solution on the patch P with Dirichlet boundary data
        %   given by BC.
        %
        %   [U, D] = SOLVE(P, BC) returns instead cell arrays containing
        %   the solution coefficients U and a vector D containing the
        %   domains of the subpatches.

            % Evaluate the solution operator for the patch:
            u = P.S * [bc ; 1]; % The 1 accounts for the particular part.

            % Construct boundary conditions for children and recurse.

            % Construct boundary indicies:
            i1 = 1:numel(P.idx1{1});
            i2 = 1:numel(P.idx2{1});
            if ( ~isempty(i1) )
                i2 = i2 + i1(end);
            end
            % CAT() is 10x faster than CELL2MAT().
            idx1 = cat(1, P.idx1{:}); % idx1 = cell2mat(p.idx1.');
            idx2 = cat(1, P.idx2{:}); % idx2 = cell2mat(p.idx2.');

            % Assemble boundary conditions for child patches:
            % - Apply the global-to-local map (i.e., the transpose of the
            %   local-to-global map) to obtain the local boundary
            %   conditions.
            ubc1 = ones(size(P.child1.S, 2)-1, 1);
            ubc1(idx1) = [bc(i1) ; P.l2g1.'*u];
            ubc2 = ones(size(P.child2.S, 2)-1, 1);
            ubc2(idx2) = [bc(i2) ; P.l2g2.'*u];

            % Solve for the child patches:
            [u1, d1] = solve(P.child1, ubc1);
            [u2, d2] = solve(P.child2, ubc2);

            % Concatenate for output:
            u = [u1 ; u2]; 
            d = [d1 ; d2];

            if ( nargout == 1 )
                % If a single output is requested, return an ultraSEMSol
                % object.
                u = ultraSEMSol(u, d);
            end

        end

        function P = updateRHS(P, rhs)
        %UPDATERHS   Update RHS of an ultraSEMParent object.
        %   P = UPDATERHS(P, RHS) replaces the existing RHS of an
        %   initialized ultraSEMParent object P with that given in RHS,
        %   which may be a constant or a function handle.

            % Update RHS of children:
            a = updateRHS(P.child1, rhs);
            b = updateRHS(P.child2, rhs);

            i1 = P.idx1{1};
            s1 = P.idx1{2};
            i2 = P.idx2{1};
            s2 = P.idx2{2};
            l2g1 = P.l2g1;
            l2g2 = P.l2g2;

            % Extract D2N maps:
            D2Na = a.D2N; D2Nb = b.D2N;
            % and discard from children
            % a.D2N = []; b.D2N = [];

            % Compute new solution operator:
            S = lsqminnorm( -( l2g1*D2Na(s1,s1)*l2g1.' + l2g2*D2Nb(s2,s2)*l2g2.' ), ...
                               l2g1*D2Na(s1,end) + l2g2*D2Nb(s2,end) );
            %                 |------------------ rhs --------------|

            % Compute new D2N maps:
            %      |--- rhs ----|
            D2N = [ D2Na(i1,end) ;
                    D2Nb(i2,end) ] ...
                + [ D2Na(i1,s1)*l2g1.' ; D2Nb(i2,s2)*l2g2.' ] * S;

            P.S(:,end) = S;
            P.D2N(:,end) = D2N;

            P.child1 = a;
            P.child2 = b;

        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STATIC METHODS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods ( Access = public, Static = true )

    end

end