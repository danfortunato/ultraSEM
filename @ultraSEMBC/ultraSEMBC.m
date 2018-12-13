classdef ultraSEMBC
%ULTRASEMBC  Boundary condition object from the ULTRASEM system.
%   An ULTRASEMBC defines the boundary conditions associated with a PDE.

% Copyright 2018 by Nick Hale and Dan Fortunato.

    %#ok<*PROP>

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties ( Access = public )

        coeffs = [] % Coefficients of Dirichlet boundary data.

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods

        function obj = ultraSEMBC(S, bc)
        %ULTRASEMBC   Class constructor for the @ultraSEMBC class.
        %   P = ultraSEMBC(domain, bc) assigns each of the inputs to their
        %   associated properties in the ultraSEMBC object P.

            % Construct empty BC:
            if ( nargin == 0 )
                return
            end

            if ( numel(S.patches) ~= 1 )
                error('ULTRASEM:ULTRASEMBC:notBuilt', ...
                    'The ultraSEM object %s must be built before creating BCs.', ...
                    inputname(1));
            end

            % Get boundary points:
            xy = S.patches{1}.xy;

            % Initialize boundary data. We should have n-2 coefficients for
            % each boundary.
            obj.coeffs = cell(size(S.patches{1}.xy));
            if ( ~isnumeric(bc) )
                % Evaluate the RHS if given a function handle:
                for k = 1:size(obj.coeffs, 1)
                    vals = feval(bc, xy{k}(:,1), xy{k}(:,2));
                    % Convert from values to coeffs:
                    obj.coeffs{k} = chebtech2.vals2coeffs(vals);
                    obj.coeffs{k} = obj.coeffs{k}(1:end-2);
                end
            elseif ( isscalar(bc) )
                % Convert a scalar to coeffs:
                for k = 1:size(obj.coeffs, 1)
                    n = size(xy{k}, 1);
                    obj.coeffs{k} = [bc ; zeros(n-3, 1)];
                end
            else
                error('ULTRASEM:ULTRASEMBC:badBC', ...
                    'Cannot evaluate boundary data.');
            end

            % CAT() is 10x faster than CELL2MAT().
            obj.coeffs = cat(1, obj.coeffs{:});

        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% METHODS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods ( Access = public, Static = false )

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STATIC METHODS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods ( Access = public, Static = true )

    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% METHODS IMPLEMENTED IN THIS FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
