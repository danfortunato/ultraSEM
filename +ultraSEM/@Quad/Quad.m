classdef Quad < ultraSEM.Mapping
%ULTRASEM.QUAD   Quadrilateral mapping object from the ULTRASEM system.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties ( Access = public )
        % Coefficients in the bilinear map:
        %   x = a1 + b1*r + c1*s + d1*r.*s
        %   y = a2 + b2*r + c2*s + d2*r.*s
        a1
        b1
        c1
        d1
        a2
        b2
        c2
        d2
        scl = @(r,s) 1+0*r; % Default scale factor for RHS
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods

        function obj = Quad( v )

            % Construct empty ULTRASEM.QUAD:
            if ( nargin == 0 )
                return
            end

            if ( ~iscell(v) )
                v = {v};
            end
            nobj = numel(v);
            obj = repelem(obj,nobj,1);
            for k = 1:nobj
                v{k} = ultraSEM.Quad.assertIsQuad(v{k});
                obj(k).v = v{k};
                if ( ~obj.isRect() )
                    obj(k).scl = @(r,s) obj(k).det(r,s).^3;
                end
            end

        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods

        function x = x(obj, r, s)
            % Map from square (r,s) to quad (x,y)
            x = obj.a1 + obj.b1*r + obj.c1*s + obj.d1*r.*s;
        end

        function y = y(obj, r, s)
            % Map from square (r,s) to quad (x,y)
            y = obj.a2 + obj.b2*r + obj.c2*s + obj.d2*r.*s;
        end

        function r = r(obj, x, y)
            % Map from quad (x,y) to square (r,s)
            A = obj.b1*obj.d2 - obj.d1*obj.b2;
            B = @(x,y) obj.b1*obj.c2 - obj.c1*obj.b2 + ...
                       obj.a1*obj.d2 - obj.d1*obj.a2 - obj.d2*x + obj.d1*y;
            C = @(x,y) obj.a1*obj.c2 - obj.c1*obj.a2 - obj.c2*x + obj.c1*y;
            if ( abs(A) < eps )
                r = -C(x,y)./B(x,y);
            else
                r = (-B(x,y) + sqrt(B(x,y).^2-4*A*C(x,y))) / (2*A);
            end
        end

        function s = s(obj, x, y)
            % Map from quad (x,y) to square (r,s)
            A = obj.c1*obj.d2 - obj.d1*obj.c2;
            B = @(x,y) obj.c1*obj.b2 - obj.b1*obj.c2 + ...
                       obj.a1*obj.d2 - obj.d1*obj.a2 - obj.d2*x + obj.d1*y;
            C = @(x,y) obj.a1*obj.b2 - obj.b1*obj.a2 - obj.b2*x + obj.b1*y;
            if ( abs(A) < eps )
                s = -C(x,y)./B(x,y);
            else
                s = (-B(x,y) - sqrt(B(x,y).^2-4*A*C(x,y))) / (2*A);
            end
        end

        function dxdr = dxdr(obj, r, s)
            dxdr = obj.b1 + obj.d1*s;
        end

        function dxds = dxds(obj, r, s)
            dxds = obj.c1 + obj.d1*r;
        end

        function dydr = dydr(obj, r, s)
            dydr = obj.b2 + obj.d2*s;
        end

        function dyds = dyds(obj, r, s)
            dyds = obj.c2 + obj.d2*r;
        end

        function drdx = drdx(obj, x, y)
            drdx = obj.dyds(x,y);
        end

        function drdy = drdy(obj, x, y)
            drdy = -obj.dxds(x,y);
        end

        function dsdx = dsdx(obj, x, y)
            dsdx = -obj.dydr(x,y);
        end

        function dsdy = dsdy(obj, x, y)
            dsdy = obj.dxdr(x,y);
        end

        function det = det(obj, r, s)
            det = (obj.b1*obj.c2 - obj.b2*obj.c1) + obj.ddetdr*r + obj.ddetds*s;
        end

        function ddetdr = ddetdr(obj)
            ddetdr = obj.b1*obj.d2 - obj.b2*obj.d1;
        end

        function ddetds = ddetds(obj)
            ddetds = obj.c2*obj.d1 - obj.c1*obj.d2;
        end

        function ddetdx = ddetdx(obj, x, y)
            ddetdx = obj.ddetdr*obj.drdx(x,y) + obj.ddetds*obj.dsdx(x,y);
        end

        function ddetdy = ddetdy(obj, x, y)
            ddetdy = obj.ddetdr*obj.drdy(x,y) + obj.ddetds*obj.dsdy(x,y);
        end

        function drdxx = drdxx(obj, x, y)
            drdxx = obj.d2*obj.drdx(x,y).*obj.det(x,y) - obj.drdx(x,y).*obj.ddetdx(x,y);
        end

        function dsdxx = dsdxx(obj, x, y)
            dsdxx = -obj.d2*obj.dsdx(x,y).*obj.det(x,y) - obj.dsdx(x,y).*obj.ddetdx(x,y);
        end

        function drdyy = drdyy(obj, x, y)
            drdyy = -obj.d1*obj.drdy(x,y).*obj.det(x,y) - obj.drdy(x,y).*obj.ddetdy(x,y);

        end

        function dsdyy = dsdyy(obj, x, y)
            dsdyy = obj.d1*obj.dsdy(x,y).*obj.det(x,y) - obj.dsdy(x,y).*obj.ddetdy(x,y);
        end

        function drdxy = drdxy(obj, x, y)
            drdxy = obj.d2*obj.drdy(x,y).*obj.det(x,y) - obj.drdx(x,y).*obj.ddetdy(x,y);
        end

        function dsdxy = dsdxy(obj, x, y)
            dsdxy = obj.d1*obj.dsdx(x,y).*obj.det(x,y) - obj.dsdy(x,y).*obj.ddetdx(x,y);
        end

        function out = dinvT11(obj, r, s)
            % Fill in Aaron's formulas for the quad transformations:
            A = obj.b1*obj.d2 - obj.d1*obj.b2;
            B = @(r,s) obj.b1*obj.c2 - obj.c1*obj.b2 + ...
                       obj.a1*obj.d2 - obj.d1*obj.a2 - obj.d2*r + obj.d1*s;
            C = @(r,s) obj.a1*obj.c2 - obj.c1*obj.a2 - obj.c2*r + obj.c1*s;
            if ( A == 0 )
                out = -(-obj.c2*(B(r,s)+obj.d2*r) + (C(r,s)+obj.c2*r)*obj.d2)./(B(r,s).^2);
            else
                out = obj.d2/(2*A) + (-2*obj.d2*B(r,s)+4*obj.c2*A)./(4*A*(sqrt(B(r,s).^2-4*A*C(r,s))));
            end
        end

        function out = dinvT12(obj, r, s)
            % Fill in Aaron's formulas for the quad transformations:
            A = obj.b1*obj.d2 - obj.d1*obj.b2;
            B = @(r,s) obj.b1*obj.c2 - obj.c1*obj.b2 + ...
                       obj.a1*obj.d2 - obj.d1*obj.a2 - obj.d2*r + obj.d1*s;
            C = @(r,s) obj.a1*obj.c2 - obj.c1*obj.a2 - obj.c2*r + obj.c1*s;
            if ( A == 0 )
                out = -(obj.c1*(B(r,s)-obj.d1*s) - (C(r,s)-obj.c1*s)*obj.d1)./(B(r,s).^2);
            else
                out = -obj.d1/(2*A) + (2*obj.d1*B(r,s)-4*obj.c1*A)./(4*A*(sqrt(B(r,s).^2-4*A*C(r,s))));
            end
        end

        function out = dinvT21(obj, r, s)
            % Fill in Aaron's formulas for the quad transformations:
            A = obj.c1*obj.d2 - obj.d1*obj.c2;
            B = @(r,s) obj.c1*obj.b2 - obj.b1*obj.c2 + ...
                       obj.a1*obj.d2 - obj.d1*obj.a2 - obj.d2*r + obj.d1*s;
            C = @(r,s) obj.a1*obj.b2 - obj.b1*obj.a2 - obj.b2*r + obj.b1*s;
            if ( A == 0 )
                out = -(-obj.b2*(B(r,s)+obj.d2*r) + (C(r,s)+obj.b2*r)*obj.d2)./(B(r,s).^2);
            else
                out = obj.d2/(2*A) - (-2*obj.d2*B(r,s)+4*obj.b2*A)./(4*A*(sqrt(B(r,s).^2-4*A*C(r,s))));
            end
        end

        function out = dinvT22(obj, r, s)
            % Fill in Aaron's formulas for the quad transformations:
            A = obj.c1*obj.d2 - obj.d1*obj.c2;
            B = @(r,s) obj.c1*obj.b2 - obj.b1*obj.c2 + ...
                       obj.a1*obj.d2 - obj.d1*obj.a2 - obj.d2*r + obj.d1*s;
            C = @(r,s) obj.a1*obj.b2 - obj.b1*obj.a2 - obj.b2*r + obj.b1*s;
            if ( A == 0 )
                out = -(obj.b1*(B(r,s)-obj.d1*s) - (C(r,s)-obj.b1*s)*obj.d1)./(B(r,s).^2);
            else
                out = -obj.d1/(2*A) - (2*obj.d1*B(r,s)-4*obj.b1*A)./(4*A*(sqrt(B(r,s).^2-4*A*C(r,s))));
            end
        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STATIC METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods ( Static )

        % Check we have valid vertices for an ULTRASEM.QUAD.
        v = assertIsQuad(v);

    end

end
