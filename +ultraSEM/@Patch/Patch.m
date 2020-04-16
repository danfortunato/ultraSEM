classdef ( Abstract ) Patch
%ULTRASEM.PATCH   Abstract patch object from the ULTRASEM system.

% Copyright 2018 by Nick Hale and Dan Fortunato.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties ( Access = public )

        domain % Domain of the patch.
        S      % Solution operator for patch.
        D2N    % Dirichlet-to-Neumann map for patch.
        edges  % Boundary edges of patch.

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% ABSTRACT METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods ( Abstract )

        % Solve a patch.
        [u, d] = solve(P, bc);

        % Update RHS of a patch.
        f = updateRHS(f, rhs);

        % Number of degrees of freedom in a patch.
        N = numel(P);

    end

end
