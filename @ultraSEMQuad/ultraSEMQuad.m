classdef ultraSEMQuad < ultraSEMMapping
%ULTRASEMQUAD   Quadrilateral mapping object from ULTRASEM system.

    %#ok<*PROP>

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        function obj = ultraSEMQuad( v )
            %ULTRASEMQUAD  Class constructor for the @ultraSEMQuad class.

            % Construct empty ultraSEMQuad:
            if ( nargin == 0 )
                return
            end
            
            % Ensure v is of the form [x, y], not [x ; y]:
            if ( size(v,1) == 2 ), v = v.'; end

            % Ensure vertices are oriented in an anticlockwise direction:
            if ( ultraSEMDomain.isClockwise(v) )
                % Switch the second and fourth indices.
                v([2,4],:) = v([4,2],:);
            end
<<<<<<< Updated upstream

=======
            
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
        
        function [Q, mergeIdx] = refine(Q, m)

            mergeIdx = {}; 
            if ( nargin == 1 ), m = 1; end
            if ( m == 0 ), return, end
            
            for j = 1:m
                nQ = numel(Q);
                Q2 = cell(nQ, 1);
                for k = 1:nQ
                    Q2{k} = refine1(Q(k));
                end
                Q = vertcat(Q2{:});
            end
            
            
        end
        
        
        function [out, mergeIdx] = refine1(Q)
            v = vertices(Q);
            c = centroid(Q).';
            vnew = (v([1 2 3 4],:) + v([2 3 4 1],:))/2;
            
            out(4,1) = ultraSEMQuad();  
            out(1) = ultraSEMQuad([v(1,:) ; vnew(1,:) ; c ; vnew(4,:)]);
            out(2) = ultraSEMQuad([v(2,:) ; vnew(2,:) ; c ; vnew(1,:)]);
            out(3) = ultraSEMQuad([v(3,:) ; vnew(3,:) ; c ; vnew(2,:)]);
            out(4) = ultraSEMQuad([v(4,:) ; vnew(4,:) ; c ; vnew(3,:)]);
            mergeIdx = {[1 2 ; 3 4], [1 2]};
            out = ultraSEMDomain(out, mergeIdx);
        end

        function Q = refineCorner(Q, k)
            if ( k ~= 1)
                error ('only k = 1 implemented')
            end
            v = vertices(Q);
            c = centroid(Q).';
            vnew = (v(1:end,:) + v([2:end, 1],:))/2;
            Q(3,1) = ultraSEMQuad();
            Q(1) = ultraSEMQuad([v(1,:) ; vnew(1,:) ; c ; vnew(4,:)]);
            Q(2) = ultraSEMQuad([v(2,:) ; v(3,:) ; c ; vnew(1,:)]);            
            Q(3) = ultraSEMQuad([v(4,:) ; vnew(4,:) ; c ; v(3,:)]);
        end
        
        function out = T1(obj, x, y)
            
            v = obj.v;
            
            % Compute the change of variables:
            M = [1 -1 -1  1;    % (-1, 1)
                 1  1 -1 -1;    % (1, -1)
                 1  1  1  1;    % (1,  1)
                 1 -1  1 -1];

            % Solve 4x4 linear system:
            c = M \ v(:,1);
            out = c(1) + c(2)*x + c(3)*y + c(4)*x.*y;
        end
        
        function out = T2(obj, x, y)
            
            v = obj.v;
            
            % Compute the change of variables:
            M = [1 -1 -1  1;    % (-1, 1)
                 1  1 -1 -1;    % (1, -1)
                 1  1  1  1;    % (1,  1)
                 1 -1  1 -1];

            % Solve 4x4 linear system:
            d = M \ v(:,2);
            out = d(1) + d(2)*x + d(3)*y + d(4)*x.*y;
        end
        
        function out = invT1(obj, s, t)
            
            v = obj.v;
            
>>>>>>> Stashed changes
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
            A = c(2)*d(4) - c(4)*d(2);
            B = @(s,t) c(2)*d(3) + c(1)*d(4) - c(4)*d(1) - c(3)*d(2) - d(4)*s + c(4)*t;
            C = @(s,t) c(1)*d(3) - c(3)*d(1) - d(3)*s + c(3)*t;
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
        
        function Q = refine(Q)
            v = vertices(Q);
            c = centroid(Q);
            vnew = (v(:,[1 2 3 4]) + v(:,[2 3 4 1]))/2;
            
            Q(4,1) = ultraSEMQuad();
            Q(1) = ultraSEMQuad([v(:,1) vnew(:,1) c vnew(:,4)]');
            Q(2) = ultraSEMQuad([v(:,2) vnew(:,2) c vnew(:,1)]');
            Q(3) = ultraSEMQuad([v(:,3) vnew(:,3) c vnew(:,2)]');
            Q(4) = ultraSEMQuad([v(:,4) vnew(:,4) c vnew(:,3)]');
        end

        function Q = refineCorner(Q, k)
            v = vertices(Q);
            c = centroid(Q);
            vnew = (v(:,1:end) + v(:,[2:end, 1]))/2;
            Q(3,1) = ultraSEMQuad();
            Q(1) = ultraSEMQuad([v(:,1) vnew(:,1) c vnew(:,4)]);
            Q(2) = ultraSEMQuad([v(:,2) v(:,3) c vnew(:,1)]);            
            Q(3) = ultraSEMQuad([v(:,4) vnew(:,4) c v(:,3)]);
        end
            
        
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STATIC METHODS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods ( Access = public, Static = true )

    end

end
