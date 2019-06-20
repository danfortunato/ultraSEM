classdef quad < ultraSEMMapping
%QUAD   Quadrilateral mapping object from ULTRASEM system.

    %#ok<*PROP>

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        function obj = quad( v )
            %QUAD  Class constructor for the @quad class.

            % Construct empty quad:
            if ( nargin == 0 )
                return
            end
            
            % Ensure v is of the form [x, y], not [x ; y]:
            if ( size(v,1) == 2 ), v = v.'; end

            % TODO: Ensure vertices are oriented in an anticlockwise
            % direction (or fix things so that clockwise is supported).

            % Compute the change of variables:
            M = [1 -1 -1  1;    % (-1, 1)
                 1  1 -1 -1;    % (1, -1)
                 1  1  1  1;    % (1,  1)
                 1 -1  1 -1];

            % Solve 4x4 linear system:
            c = M \ v(:,1);
            obj.T1 = @(x,y) c(1) + c(2)*x + c(3)*y + c(4)*x.*y;
            d = M \ v(:,2);
            obj.T2 = @(x,y) d(1) + d(2)*x + d(3)*y + d(4)*x.*y;

            % Fill in Aaron's formulas for the quad transformations:
            A = c(2)*d(4)-c(4)*d(2);
            B = @(s,t) c(2)*d(3)+c(1)*d(4)-c(4)*d(1)-c(3)*d(2)-d(4)*s+c(4)*t;
            C = @(s,t) c(1)*d(3)-c(3)*d(1)-d(3)*s+c(3)*t;
            if ( A == 0 )
                obj.invT1 = @(s,t) -C(s,t)./B(s,t);
                obj.dinvT1{1} = @(s,t) -(-d(3)*(B(s,t)+d(4)*s) + (C(s,t)+d(3)*s)*d(4))./(B(s,t).^2);
                obj.dinvT1{2} = @(s,t) -(c(3)*(B(s,t)-c(4)*t) - (C(s,t)-c(3)*t)*c(4))./(B(s,t).^2);
                obj.d2invT1{1} = @(s,t) (2*d(4)*(-c(2)*d(3)^2 + c(4)*d(3)*(d(1) - t) + c(3)*(d(2)*d(3) - d(1)*d(4) + d(4)*t)))./(c(3)*d(2) - c(2)*d(3) - c(1)*d(4) + d(4)*s + c(4)*(d(1) - t)).^3;
                obj.d2invT1{2} = @(s,t) (c(3)*d(4)*(-c(3)*d(2) + c(2)*d(3) + c(1)*d(4) - d(4)*s) + c(4)^2*d(3)*(-d(1) + t) + c(4)*(d(3)*(c(2)*d(3) - c(1)*d(4) + d(4)*s) - c(3)*(d(2)*d(3) - d(1)*d(4) + d(4)*t)))./(c(3)*d(2) - c(2)*d(3) - c(1)*d(4) + d(4)*s + c(4)*(d(1) - t)).^3;
                obj.d2invT1{3} = @(s,t) (2*c(4)*(c(3)^2*d(2) - c(3)*(c(2)*d(3) + d(4)*(c(1) - s)) + c(4)*d(3)*(c(1) - s)))./(c(3)*d(2) - c(2)*d(3) - c(1)*d(4) + d(4)*s + c(4)*(d(1) - t)).^3;
            else
                obj.invT1 = @(s,t) (-B(s,t) + sqrt(B(s,t).^2-4*A*C(s,t)))/(2*A);
                obj.dinvT1{1} = @(s,t) d(4)/(2*A) + (-2*d(4)*B(s,t)+4*d(3)*A)./(4*A*(sqrt(B(s,t).^2-4*A*C(s,t))));
                obj.dinvT1{2} = @(s,t) -c(4)/(2*A) + (2*c(4)*B(s,t)-4*c(3)*A)./(4*A*(sqrt(B(s,t).^2-4*A*C(s,t))));
                obj.d2invT1{1} = @(s,t) (2*(c(4)*d(3) - c(3)*d(4))*(d(2)*d(3) + d(4)*(-d(1) + t)))./((c(4)*d(1) + c(3)*d(2)-c(2)*d(3) - c(1)*d(4) + d(4)*s - c(4)*t).^2 + 4*(c(4)*d(2) - c(2)*d(4))*(d(3)*(c(1) - s) + c(3)*(-d(1) + t))).^(3/2);
                obj.d2invT1{2} = @(s,t) ((-c(4)*d(3) + c(3)*d(4))*(c(3)*d(2) + c(2)*d(3) - c(1)*d(4) + d(4)*s + c(4)*(-d(1) + t)))./((c(4)*d(1) + c(3)*d(2) - c(2)*d(3) - c(1)*d(4) + d(4)*s - c(4)*t).^2 + 4*(c(4)*d(2) - c(2)*d(4)).*(d(3)*(c(1) - s) + c(3)*(-d(1) + t))).^(3/2);
                obj.d2invT1{3} = @(s,t) (2*(c(4)*d(3) - c(3)*d(4))*(c(2)*c(3) + c(4)*(-c(1) + s)))./((c(4)*d(1) + c(3)*d(2) - c(2)*d(3) - c(1)*d(4) + d(4)*s - c(4)*t).^2 + 4*(c(4)*d(2) - c(2)*d(4)).*(d(3)*(c(1) - s) + c(3)*(-d(1) + t))).^(3/2);
            end

            A = c(3)*d(4)-c(4)*d(3);
            B = @(s,t) c(3)*d(2)+c(1)*d(4)-c(4)*d(1)-c(2)*d(3)-d(4)*s+c(4)*t;
            C = @(s,t) c(1)*d(2)-c(2)*d(1)-d(2)*s+c(2)*t;
            if ( A == 0 )
                obj.invT2 = @(s,t) -C(s,t)./B(s,t);
                obj.dinvT2{1} = @(s,t) -(-d(2)*(B(s,t)+d(4)*s) + (C(s,t)+d(2)*s)*d(4))./(B(s,t).^2);
                obj.dinvT2{2} = @(s,t) -(c(2)*(B(s,t)-c(4)*t) - (C(s,t)-c(2)*t)*c(4))./(B(s,t).^2);
                obj.d2invT2{1} = @(s,t) (2*d(4)*(-c(3)*d(2)^2 + c(4)*d(2)*(d(1) - t) + c(2)*(d(3)*d(2) - d(1)*d(4) + d(4)*t)))./(c(2)*d(3) - c(3)*d(2) - c(1)*d(4) + d(4)*s + c(4)*(d(1) - t)).^3;
                obj.d2invT2{2} = @(s,t) (c(2)*d(4)*(-c(2)*d(3) + c(3)*d(2) + c(1)*d(4) - d(4)*s) + c(4)^2*d(2)*(-d(1) + t) + c(4)*(d(2)*(c(3)*d(2) - c(1)*d(4) + d(4)*s) - c(2)*(d(3)*d(2) - d(1)*d(4) + d(4)*t)))./(c(2)*d(3) - c(3)*d(2) - c(1)*d(4) + d(4)*s + c(4)*(d(1) - t)).^3;
                obj.d2invT2{3} = @(s,t) (2*c(4)*(c(2)^2*d(3) - c(2)*(c(3)*d(2) + d(4)*(c(1) - s)) + c(4)*d(2)*(c(1) - s)))./(c(2)*d(3) - c(3)*d(2) - c(1)*d(4) + d(4)*s + c(4)*(d(1) - t)).^3;
            else
                obj.invT2 = @(s,t) (-B(s,t) - sqrt(B(s,t).^2-4*A*C(s,t)))/(2*A);
                obj.dinvT2{1} = @(s,t) d(4)/(2*A) - (-2*d(4)*B(s,t)+4*d(2)*A)./(4*A*(sqrt(B(s,t).^2-4*A*C(s,t))));
                obj.dinvT2{2} = @(s,t) -c(4)/(2*A) - (2*c(4)*B(s,t)-4*c(2)*A)./(4*A*(sqrt(B(s,t).^2-4*A*C(s,t))));
                obj.d2invT2{1} = @(s,t) (2*(-c(4)*d(2) + c(2)*d(4))*(d(3)*d(2) + d(4)*(-d(1) + t)))./((c(4)*d(1) + c(2)*d(3)-c(3)*d(2) - c(1)*d(4) + d(4)*s - c(4)*t).^2 + 4*(c(4)*d(3) - c(3)*d(4))*(d(2)*(c(1) - s) + c(2)*(-d(1) + t))).^(3/2);
                obj.d2invT2{2} = @(s,t) -((-c(4)*d(2) + c(2)*d(4))*(c(2)*d(3) + c(3)*d(2) - c(1)*d(4) + d(4)*s + c(4)*(-d(1) + t)))./((c(4)*d(1) + c(2)*d(3) - c(3)*d(2) - c(1)*d(4) + d(4)*s - c(4)*t).^2 + 4*(c(4)*d(3) - c(3)*d(4)).*(d(2)*(c(1) - s) + c(2)*(-d(1) + t))).^(3/2);
                obj.d2invT2{3} = @(s,t) (2*(-c(4)*d(2) + c(2)*d(4))*(c(3)*c(2) + c(4)*(-c(1) + s)))./((c(4)*d(1) + c(2)*d(3) - c(3)*d(2) - c(1)*d(4) + d(4)*s - c(4)*t).^2 + 4*(c(4)*d(3) - c(3)*d(4)).*(d(2)*(c(1) - s) + c(2)*(-d(1) + t))).^(3/2);
            end

        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% METHODS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods
        
        function c = centroid(Q)
            v = vertices(Q);
            [xmid, ymid] = centroid(polyshape(v.'));
            c = [xmid ; ymid]; 
        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STATIC METHODS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods ( Access = public, Static = true )

    end

end
