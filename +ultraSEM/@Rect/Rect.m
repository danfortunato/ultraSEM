classdef Rect < ultraSEM.Quad
%ULTRASEM.RECT   Rectangular mapping object from the ULTRASEM system.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods

        function obj = Rect( dom )

            % Construct empty ULTRASEM.RECT:
            args = {};
            if ( nargin ~= 0 )
                args{1} = util.rect2quad(dom);
            end

            % Call the superclass constructor
            obj = obj@ultraSEM.Quad(args{:});

        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STATIC METHODS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods ( Access = private, Static = true )

        % Check we have valid vertices for an ULTRASEM.RECT.
        v = assertIsRect( v );

    end

end
