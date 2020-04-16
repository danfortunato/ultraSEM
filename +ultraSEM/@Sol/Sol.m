classdef Sol
%ULTRASEM.SOL   Solution object from the ULTRASEM system.
%   An ULTRASEM.SOL is the result of solving a PDE using ULTRASEM. It
%   allows elementary exploration of the solution, such as various methods 
%   of plotting and evaluation within the domain.

% Copyright 2018 by Nick Hale and Dan Fortunato.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties ( Access = public )

        domain % Domain of the patch.
        u      % Coefficients of solution.

    end

    properties ( Constant )

        nplotpts = 200; % Number of plot points.

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods

        function obj = Sol(u, d)

            % Construct empty patch:
            if ( nargin == 0 )
                return
            end

            obj.domain = d;
            obj.u = u;

        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods

        function n = numArgumentsFromSubscript(obj,s,indexingContext)
        %NUMARGUMENTSFROMSUBSCRIPT   Number of arguments for customized indexing methods.
        %   Overloading NUMEL() gives the wrong NARGOUT for SUBSREF().
        %   Defining this function fixes it.
        %
        % See also NUMEL, NARGOUT, SUBSREF.

            n = 1;

        end

    end

end
