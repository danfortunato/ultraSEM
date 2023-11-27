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

        function obj = BC(val, Dir, neu)

            if ( nargin == 0 )
                % Return an empty object:
                return
            end

            obj.dir = Dir;
            obj.neu = neu;
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
        %   VAL. VAL may be a constant or a function handle of two
        %   variables (i.e. VAL = @(x,y) ...).
        %
        %   ULTRASEM.BC.DIRICHLET() assumes VAL = 0.
        
            % Assume zero if not specified:
            if ( nargin == 0 ), val = 0; end
            
            bc = ultraSEM.BC(val, 1, 0);
        end

        function bc = neumann(val)
        %NEUMANN   Construct Neumann boundary conditions.
        %   ULTRASEM.BC.NEUMANN(VAL) constructs an ULTRASEM.BC representing
        %   a Neumann boundary condition with boundary data VAL. VAL may be
        %   a constant or a function handle of two variables (i.e., VAL =
        %   @(x,y) ...).
        %
        %   ULTRASEM.BC.NEUMANN() assumes VAL = 0.
        
            % Assume zero if not specified:
            if ( nargin == 0 ), val = 0; end
        
            bc = ultraSEM.BC(val, 0, 1);
        end

        function bc = robin(val, a, b)
        %ROBIN   Construct Robin boundary conditions.
        %   ULTRASEM.BC.ROBIN(VAL, A, B) constructs an ULTRASEM.BC
        %   representing the Robin boundary condition A + B*d/dn with
        %   boundary data VAL and constants A and B. VAL, A, and B may be a
        %   constants or a function handles of two variables (i.e., VAL =
        %   @(x,y) ...).
        
            bc = ultraSEM.BC(val, a, b);
        end
        
        function bc = getBCsFromPlot(h)
            bc = ultraSEM.BC;
            for k = 1:numel(h)
                data = get(h(k), 'UserData');
                data2 = cellfun(@str2double, data, 'uniformoutput', false);
                if ( isnan(data2{3}) )
                    data2{3} = str2func(['@(x,y) ' vectorize(data{3})]);
                end
                bc(k) = ultraSEM.BC(data2{3}, data2{1}, data2{2});
            end
        end

    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods
        
        function out = validBC(bc)
            %VALIDBC  Check a BC is contains all the necessary information.
            
            if ( numel(bc) > 1 )
                out = zeros(size(bc));
                for k = 1:numel(bc)
                    out(k) = validBC(bc(k));
                end
                return
            end
            
            out =  ( ~isempty(bc.dir) && ~isempty(bc.neu) && ...
                ~isempty(bc.val) && (bc.dir ~= 0 || bc.neu ~= 0) );
            
        end
        
    end

    

    

end
