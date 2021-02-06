classdef BC < matlab.mixin.Heterogeneous
%ULTRASEM.BC   Boundary condition object from the ULTRASEM system.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % The properties in the BC object correspond to the boundary condition:
    %
    %    dir(k)*u + neu(k)*du/dn = val

    properties ( Access = public )

        dir
        neu
        val

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods

        function obj = BC(val, scl)

            if ( nargin == 0 )
                return
            end

            if ( nargin == 1 )
                % Default to Dirichlet boundary conditions:
                scl = [1 0];
            end

            if ( isnumeric(scl) && (isscalar(scl) || (isvector(scl) && length(scl)==2) ) )
                if ( isscalar(scl) )
                    obj.dir = scl;
                    obj.neu = scl;
                else
                    obj.dir = scl(1);
                    obj.neu = scl(2);
                end
            else
                error('ULTRASEM:BC:bc:invalid', ...
                    'Invalid constants in boundary condition.');
            end

            obj.val = val;

        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STATIC METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods ( Access = public, Static = true )

        function bc = dirichlet(val)
        %DIRICHLET   Construct Dirichlet boundary conditions.
        %   ULTRASEM.BC.DIRICHLET(VAL) constructs an ULTRASEM.BC
        %   representing a Dirichlet boundary condition with boundary data
        %   VAL.
            bc = ultraSEM.BC(val, [1 0]);
        end

        function bc = neumann(val)
        %NEUMANN   Construct Neumann boundary conditions.
        %   ULTRASEM.BC.NEUMANN(VAL) constructs an ULTRASEM.BC representing
        %   a Neumann boundary condition with boundary data VAL.
            bc = ultraSEM.BC(val, [0 1]);
        end

        function bc = robin(val, a, b)
        %ROBIN   Construct Robin boundary conditions.
        %   ULTRASEM.BC.ROBIN(VAL, A, B) constructs an ULTRASEM.BC
        %   representing the Robin boundary condition A + B*d/dn with
        %   boundary data VAL and constants A and B.
            bc = ultraSEM.BC(val, [a b]);
        end

    end

end
