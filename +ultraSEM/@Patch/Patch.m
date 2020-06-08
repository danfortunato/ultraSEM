classdef ( Abstract ) Patch
%ULTRASEM.PATCH   Abstract patch object from the ULTRASEM system.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties ( Access = public )

        domain  % Domain of the patch.
        S       % Solution operator for patch.
        D2N     % Dirichlet-to-Neumann map for patch.
        edges   % Boundary edges of patch.
        D2N_scl % Cell array of scalars or function handles.
                % The k-th entry is the scaling for the D2N map on side k.

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

        % Number of patches in a patch.
        N = length(P);

    end

end
