classdef Tri < Mapping
%ULTRASEM.TRI   Triangle mapping object from the ULTRASEM system.
%   For now this just maps [-1,1]^2 -> right-angled triangle.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods

        function obj = Tri( v )

            % Construct empty ULTRASEM.TRI:
            if ( nargin == 0 )
                return
            end

            % Construct multiple ULTRASEM.TRIs:
            if ( iscell(v) )
                obj(numel(v),1) = ultraSEM.Tri();
                for k = 1:numel(v)
                    obj(k,1) = ultraSEM.Tri( v{k} );
                end
                return
            end

            v = ultraSEM.Tri.assertIsTri(v);

            obj.v = v;

        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% METHODS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods

        function obj = parametrize(obj)

        end

        function out = T1(obj, x, y)
            out = (x-1)/2;
        end

        function out = T2(obj, x, y)
            out = (x-1)/2.*(y-1)/2;
        end

        function out = invT1(obj, s, t)
            out = 2*s+1;
        end

        function out = invT2(obj, s, t)
            out = 2*(t./s)+1;
        end

        function out = dinvT11(obj, s, t)
            out = 2+0*s;
        end

        function out = dinvT12(obj, s, t)
            out = 0*s;
        end

        function out = dinvT21(obj, s, t)
            out = -2*(t./s.^2);
        end

        function out = dinvT22(obj, s, t)
            out = 2./s;
        end

        function out = d2invT11(obj, s, t)
            out = 0*s;
        end

        function out = d2invT12(obj, s, t)
            out = 0*s;
        end

        function out = d2invT13(obj, s, t)
            out = 0*s;
        end

        function out = d2invT21(obj, s, t)
            out = 4*(t./s.^3);
        end

        function out = d2invT22(obj, s, t)
            out = -2./s.^2;
        end

        function out = d2invT23(obj, s, t)
            out = 0*s;
        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STATIC METHODS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods ( Access = public, Static = true )

        % Check we have valid vertices for an ULTRASEM.TRI.
        v = assertIsTri(v);

    end

end
