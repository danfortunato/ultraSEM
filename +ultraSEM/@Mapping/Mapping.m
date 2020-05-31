classdef Mapping < matlab.mixin.Heterogeneous
%ULTRASEM.MAPPING   Abstract mapping object from the ULTRASEM system.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties ( Access = public )

        v

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods

        % Compute the centroid of an ULTRASEM.MAPPING.
        c = centroid(Q);

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% ABSTRACT METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods ( Abstract )

        % The parametrization defining the mapping.
        obj = parametrize(obj);

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% SEALED METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods ( Sealed )

        % Plot the mapped domain and grid.
        varargout = plot(T, varargin);

        % Plot an ULTRASEM.MAPPING as a mesh.
        varargout = mesh(T, varargin);

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %% METHODS IMPLEMENTED IN THIS FILE:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
    methods

        function obj = set.v(obj, value)
        %MAPPING.SET.V   Set vertices of a mapping.
        %   MAPPING.SET.V automatically recomputes the parametrization for
        %   the new vertices V.

            obj.v = value;
            obj = parametrize(obj);

        end

    end

    methods (Static, Sealed, Access = protected)

        function defaultObject = getDefaultScalarElement
        %GETDEFAULTSCALARELEMENT   Create a default object for heterogeneous arrays.

            defaultObject = ultraSEM.Quad;

        end

    end

end
