classdef ultraSEMParent < ultraSEMPatch

% Copyright 2018 by Nick Hale and Dan Fortunato.

    %#ok<*PROP>

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties ( Access = public )

        child1 = []     % Child patch
        child2 = []     % Child patch
        idx1            % How p.xy relates to p.child1.xy
        idx2            % How p.xy relates to p.child2.xy
        flip1           % How to flip the interface for child1
        flip2           % How to flip the interface for child2

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods

        function P = ultraSEMParent(domain, S, D2N, xy, child1, child2, idx1, idx2, flip1, flip2)
        %ULTRASEMPARENT   Class constructor for the @ultraSEMParent class.
        %   P = ultraSEMParent(DOMAIN, S, D2N, XY, CHILD1, CHILD2, IDX1, IDX2)
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
            P.xy = xy;

            % Assign children:
            P.child1 = child1;
            P.child2 = child2;
            P.idx1 = idx1;
            P.idx2 = idx2;
            P.flip1 = flip1;
            P.flip2 = flip2;

        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% METHODS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods ( Access = public, Static = false )

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
            ubc1 = ones(size(P.child1.S, 2)-1, 1);
            ubc1(idx1) = [bc(i1) ; P.flip1.*u];
            ubc2 = ones(size(P.child2.S, 2)-1, 1);
            ubc2(idx2) = [bc(i2) ; P.flip2.*u];

            % Solve for the child patches:
            [u1, d1] = solve(P.child1, ubc1);
            [u2, d2] = solve(P.child2, ubc2);

            % Concatenate for output:
            u = [u1 ; u2]; 
            
            % Cast rectangular domains to Quads.
            % TODO: This is annoying. It will be better when we don't cheat and
            % do rectangular domains properly.
            if ( isnumeric(d1) && ~isnumeric(d2) )
                for k = 1:size(d1, 1)
                    d1k = d1(k,:);
                    v1 = d1k([1 3 ; 2 3 ; 2 4 ; 1 4]);
                    q1(k,1) = ultraSEMQuad(v1);
                end
                d1 = q1;
            elseif ( ~isnumeric(d1) && isnumeric(d2) )
                for k = 1:size(d2, 1)
                    d2k = d2(k,:);
                    v2 = d2k([1 3 ; 2 3 ; 2 4 ; 1 4]);
                    q2(k,1) = ultraSEMQuad(v2);
                end
                d2 = q2;
            end
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
            flip1 = P.flip1;
            flip2 = P.flip2;

            % Extract D2N maps:
            D2Na = a.D2N; D2Nb = b.D2N;
            % and discard from children
            % a.D2N = []; b.D2N = [];

            % Compute new solution operator:
            S = lsqminnorm( -( flip1.*D2Na(s1,s1).*flip1.' + flip2.*D2Nb(s2,s2).*flip2.' ), ...
                               flip1.*D2Na(s1,end) + flip2.*D2Nb(s2,end) );
            %                 |------------------ rhs ------------------|

            % Compute new D2N maps:
            %      |--- rhs ----|
            D2N = [ D2Na(i1,end) ;
                    D2Nb(i2,end) ] ...
                + [ D2Na(i1,s1).*flip1.' ; D2Nb(i2,s2).*flip2.' ] * S;

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
