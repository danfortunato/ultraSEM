classdef ultraSEMTri < ultraSEMMapping
%ULTRASEMTRI   Triangle mapping object from ULTRASEM system.
% For now this just maps [-1,1]^2 -> right-angled triangle.

    %#ok<*PROP>

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods

        function obj = ultraSEMTri( v )
            %ULTRASEMTRI  Class constructor for the @ultraSEMTri class.

            % Construct empty ultraSEMTri:
            if ( nargin == 0 )
                return
            end

            % Construct multiple ultraSEMTris:
            if ( iscell(v) )
                obj(numel(v),1) = ultraSEMTri();
                for k = 1:numel(v)
                    obj(k,1) = ultraSEMTri( v{k} );
                end
                return
            end

            % v = ultraSEMQuad.assertIsTri(v);

            obj.v = v;

        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% METHODS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods

        function c = centroid(Q)
            [xmid, ymid] = centroid(polyshape(Q.v));
            c = [xmid ; ymid];
        end

        function out = isRect(M)
            out = false;
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

        function v = assertIsTri(v)
        %ASSERTISTRI   Check we have valid vertices for an ultraSEMTri.

            if ( ~isnumeric(v) )
                error('Input should be numeric.');
            end
            % Ensure v is of the form [x, y], not [x ; y]:
            if ( size(v,2) ~= 2 ), v = v.'; end
            % Check dimension:
            if ( size(v,2) ~= 2 || size(v, 1) ~= 3 )
                error('Incorrect vertices dimension.')
            end
            % Ensure vertices are oriented in an anticlockwise direction:
            if ( ultraSEMDomain.isClockwise(v) )
                % Switch the second and third indices.
                v([2,3],:) = v([3,2],:);
            end
        end

    end

end
