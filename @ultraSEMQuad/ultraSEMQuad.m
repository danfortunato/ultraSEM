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
            
            % Construct multiple ultraSEMQuads:
            if ( iscell(v) )
                obj(numel(v),1) = ultraSEMQuad();
                for k = 1:numel(v)
                    obj(k,1) = ultraSEMQuad( v{k} );
                end
                return
            end
            
            v = ultraSEMQuad.assertIsQuad(v);
                        
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
       
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% METHODS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods
                
        function [Q, mergeIdx] = refine(Q, m)
            
            % Parse inputs:
            mergeIdx = {}; 
            if ( nargin == 1 ), m = 1; end
            if ( m == 0 ), return, end

            n = numel(Q);       % Number of Quads
            
            % Initialize:
            v = cell(4*n, 1);
            idx1 = cell(n,1);
            idx2 = cell(n,1);
            
            % Divide into four using centroid and centre of edges:
            for k = 1:n
                vk = Q(k).v;
                c = centroid(Q(k)).';
                vnew = (vk([1 2 3 4],:) + vk([2 3 4 1],:))/2;
                v(4*(k-1)+(1:4)) = ...
                     {[vk(1,:) ; vnew(1,:) ; c ; vnew(4,:)];
                      [vnew(1,:) ; vk(2,:) ; vnew(2,:) ; c];
                      [c ; vnew(2,:) ; vk(3,:) ; vnew(3,:) ; ];
                      [vnew(4,:) ; c ; vnew(3,:) ; vk(4,:)]};                  
                idx1{k} = [1 2 ; 3 4] + 4*(k-1);
                idx2{k} = [1 2] + 2*(k-1);
            end
            
            % Assemble Quad:
            mergeIdx = {vertcat(idx1{:}), vertcat(idx2{:})};
            Q = ultraSEMQuad(v);           
            
            % Recurse for more refinement:
            [Q, idxNew] = refine(Q, m-1);
            
            % Append merge info:
            mergeIdx = [idxNew, mergeIdx];
        end        
% 
%         function Q = refineCorner(Q, k)
%             if ( k ~= 1)
%                 error ('only k = 1 implemented')
%             end
%             v = vertices(Q);
%             c = centroid(Q).';
%             vnew = (v(1:end,:) + v([2:end, 1],:))/2;
%             Q(3,1) = ultraSEMQuad();
%             Q(1) = ultraSEMQuad([v(1,:) ; vnew(1,:) ; c ; vnew(4,:)]);
%             Q(2) = ultraSEMQuad([v(2,:) ; v(3,:) ; c ; vnew(1,:)]);            
%             Q(3) = ultraSEMQuad([v(4,:) ; vnew(4,:) ; c ; v(3,:)]);
%         end

        function [Q, idx] = refinePoint(Q, z)
            loc = find(ismember(Q.v, z, 'rows'));
            if ( any(loc) )
                [Q, idx] = refineCorner(Q, loc);
            else
                idx = {[1, NaN]};
            end
        end
        
        function [Q, idx] = refineCorner(Q, k)
            % TODO: Transpose everythng.
            v = vertices(Q)';
            c = centroid(Q);
            Q(3,1) = ultraSEMQuad();
            if ( k == 1 )
                vnew = (v(:,1) + v(:,[2,4]))/2;
                v1 = [v(:,1) vnew(:,1) c vnew(:,2)];
                v2 = [v(:,[2 3]) c vnew(:,1)];
                v3 = [v(:,4) vnew(:,2) c v(:,3)];
            elseif ( k == 2 )
                vnew = (v(:,2) + v(:,[3,1]))/2;
                v1 = [v(:,2) vnew(:,1) c vnew(:,2)];
                v2 = [v(:,[3 4]) c vnew(:,1)];
                v3 = [v(:,1) vnew(:,2) c v(:,4)];                
            elseif ( k == 3 )
                vnew = (v(:,3) + v(:,[4,2]))/2;
                v1 = [v(:,3) vnew(:,1) c vnew(:,2)];
                v2 = [v(:,[4 1]) c vnew(:,1)];
                v3 = [v(:,2) vnew(:,2) c v(:,1)]; 
            elseif ( k == 4 )
                vnew = (v(:,4) + v(:,[1,3]))/2;
                v1 = [v(:,4) vnew(:,1) c vnew(:,2)];
                v2 = [v(:,[1 2]) c vnew(:,1)];
                v3 = [v(:,3) vnew(:,2) c v(:,2)];                 
            end 
            Q(3) = ultraSEMQuad(v1');
            Q(2) = ultraSEMQuad(v2');            
            Q(1) = ultraSEMQuad(v3');
            idx = {[1 2 ; 3 NaN], [1 2]};
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
            
            % Compute the change of variables:
            M = [1 -1 -1  1;    % (-1, 1)
                 1  1 -1 -1;    % (1, -1)
                 1  1  1  1;    % (1,  1)
                 1 -1  1 -1];

            % Solve 4x4 linear system:
            c = M \ v(:,1);
            d = M \ v(:,2);

            % Fill in Aaron's formulas for the quad transformations:
            A = c(2)*d(4) - c(4)*d(2);
            B = @(s,t) c(2)*d(3) + c(1)*d(4) - c(4)*d(1) - c(3)*d(2) - d(4)*s + c(4)*t;
            C = @(s,t) c(1)*d(3) - c(3)*d(1) - d(3)*s + c(3)*t;
            if ( A == 0 )
                out = -C(s,t)./B(s,t);
            else
                out = (-B(s,t) + sqrt(B(s,t).^2-4*A*C(s,t))) / (2*A);
            end
        end
        
        function out = invT2(obj, s, t)
            
            v = obj.v;
            
            % Compute the change of variables:
            M = [1 -1 -1  1;    % (-1, 1)
                 1  1 -1 -1;    % (1, -1)
                 1  1  1  1;    % (1,  1)
                 1 -1  1 -1];

            % Solve 4x4 linear system:
            c = M \ v(:,1);
            d = M \ v(:,2);

            % Fill in Aaron's formulas for the quad transformations:
            A = c(3)*d(4)-c(4)*d(3);
            B = @(s,t) c(3)*d(2) + c(1)*d(4) - c(4)*d(1) - c(2)*d(3) - d(4)*s+c(4)*t;
            C = @(s,t) c(1)*d(2) - c(2)*d(1) - d(2)*s + c(2)*t;
            if ( A == 0 )
                out = -C(s,t)./B(s,t);
            else
                out = (-B(s,t) - sqrt(B(s,t).^2-4*A*C(s,t))) / (2*A);
            end
        end
        
        function out = dinvT11(obj, s, t)
            
            v = obj.v;
            
            % Compute the change of variables:
            M = [1 -1 -1  1;    % (-1, 1)
                 1  1 -1 -1;    % (1, -1)
                 1  1  1  1;    % (1,  1)
                 1 -1  1 -1];

            % Solve 4x4 linear system:
            c = M \ v(:,1);
            d = M \ v(:,2);

            % Fill in Aaron's formulas for the quad transformations:
            A = c(2)*d(4) - c(4)*d(2);
            B = @(s,t) c(2)*d(3) + c(1)*d(4) - c(4)*d(1) - c(3)*d(2) - d(4)*s + c(4)*t;
            C = @(s,t) c(1)*d(3) - c(3)*d(1) - d(3)*s + c(3)*t;
            if ( A == 0 )
                out = -(-d(3)*(B(s,t)+d(4)*s) + (C(s,t)+d(3)*s)*d(4))./(B(s,t).^2);    
            else
                out = d(4)/(2*A) + (-2*d(4)*B(s,t)+4*d(3)*A)./(4*A*(sqrt(B(s,t).^2-4*A*C(s,t))));
            end
        end
        
        function out = dinvT12(obj, s, t)
            
            v = obj.v;
            
            % Compute the change of variables:
            M = [1 -1 -1  1;    % (-1, 1)
                 1  1 -1 -1;    % (1, -1)
                 1  1  1  1;    % (1,  1)
                 1 -1  1 -1];

            % Solve 4x4 linear system:
            c = M \ v(:,1);
            d = M \ v(:,2);

            % Fill in Aaron's formulas for the quad transformations:
            A = c(2)*d(4) - c(4)*d(2);
            B = @(s,t) c(2)*d(3) + c(1)*d(4) - c(4)*d(1) - c(3)*d(2) - d(4)*s + c(4)*t;
            C = @(s,t) c(1)*d(3) - c(3)*d(1) - d(3)*s + c(3)*t;                
            if ( A == 0 )
                out = -(c(3)*(B(s,t)-c(4)*t) - (C(s,t)-c(3)*t)*c(4))./(B(s,t).^2);
            else
                out = -c(4)/(2*A) + (2*c(4)*B(s,t)-4*c(3)*A)./(4*A*(sqrt(B(s,t).^2-4*A*C(s,t))));
            end
        end   
        
        function out = dinvT21(obj, s, t)
            
            v = obj.v;
            
            % Compute the change of variables:
            M = [1 -1 -1  1;    % (-1, 1)
                 1  1 -1 -1;    % (1, -1)
                 1  1  1  1;    % (1,  1)
                 1 -1  1 -1];

            % Solve 4x4 linear system:
            c = M \ v(:,1);
            d = M \ v(:,2);

            % Fill in Aaron's formulas for the quad transformations:
            A = c(3)*d(4)-c(4)*d(3);
            B = @(s,t) c(3)*d(2)+c(1)*d(4)-c(4)*d(1)-c(2)*d(3)-d(4)*s+c(4)*t;
            C = @(s,t) c(1)*d(2)-c(2)*d(1)-d(2)*s+c(2)*t;
            if ( A == 0 )
                out = -(-d(2)*(B(s,t)+d(4)*s) + (C(s,t)+d(2)*s)*d(4))./(B(s,t).^2);
            else
                out = d(4)/(2*A) - (-2*d(4)*B(s,t)+4*d(2)*A)./(4*A*(sqrt(B(s,t).^2-4*A*C(s,t))));
            end
        end

        function out = dinvT22(obj, s, t)
            
            v = obj.v;
            
            % Compute the change of variables:
            M = [1 -1 -1  1;    % (-1, 1)
                 1  1 -1 -1;    % (1, -1)
                 1  1  1  1;    % (1,  1)
                 1 -1  1 -1];

            % Solve 4x4 linear system:
            c = M \ v(:,1);
            d = M \ v(:,2);

            % Fill in Aaron's formulas for the quad transformations:
            A = c(3)*d(4)-c(4)*d(3);
            B = @(s,t) c(3)*d(2)+c(1)*d(4)-c(4)*d(1)-c(2)*d(3)-d(4)*s+c(4)*t;
            C = @(s,t) c(1)*d(2)-c(2)*d(1)-d(2)*s+c(2)*t;            
            if ( A == 0 )
                out = -(c(2)*(B(s,t)-c(4)*t) - (C(s,t)-c(2)*t)*c(4))./(B(s,t).^2);
            else
                out = -c(4)/(2*A) - (2*c(4)*B(s,t)-4*c(2)*A)./(4*A*(sqrt(B(s,t).^2-4*A*C(s,t))));
            end
        end   
        
        function out = d2invT11(obj, s, t)
            
            v = obj.v;
            
            % Compute the change of variables:
            M = [1 -1 -1  1;    % (-1, 1)
                 1  1 -1 -1;    % (1, -1)
                 1  1  1  1;    % (1,  1)
                 1 -1  1 -1];

            % Solve 4x4 linear system:
            c = M \ v(:,1);
            d = M \ v(:,2);

            % Fill in Aaron's formulas for the quad transformations:
            A = c(2)*d(4) - c(4)*d(2);
            if ( A == 0 )
                out = (2*d(4)*(-c(2)*d(3)^2 + c(4)*d(3)*(d(1) - t) + c(3)*(d(2)*d(3) - d(1)*d(4) + d(4)*t)))./(c(3)*d(2) - c(2)*d(3) - c(1)*d(4) + d(4)*s + c(4)*(d(1) - t)).^3;  
            else
                out = (2*(c(4)*d(3) - c(3)*d(4))*(d(2)*d(3) + d(4)*(-d(1) + t)))./((c(4)*d(1) + c(3)*d(2)-c(2)*d(3) - c(1)*d(4) + d(4)*s - c(4)*t).^2 + 4*(c(4)*d(2) - c(2)*d(4))*(d(3)*(c(1) - s) + c(3)*(-d(1) + t))).^(3/2);
            end
        end
        
        function out = d2invT12(obj, s, t)
            
            v = obj.v;
            
            % Compute the change of variables:
            M = [1 -1 -1  1;    % (-1, 1)
                 1  1 -1 -1;    % (1, -1)
                 1  1  1  1;    % (1,  1)
                 1 -1  1 -1];

            % Solve 4x4 linear system:
            c = M \ v(:,1);
            d = M \ v(:,2);

            % Fill in Aaron's formulas for the quad transformations:
            A = c(2)*d(4) - c(4)*d(2);
            if ( A == 0 )
                out = (c(3)*d(4)*(-c(3)*d(2) + c(2)*d(3) + c(1)*d(4) - d(4)*s) + c(4)^2*d(3)*(-d(1) + t) + c(4)*(d(3)*(c(2)*d(3) - c(1)*d(4) + d(4)*s) - c(3)*(d(2)*d(3) - d(1)*d(4) + d(4)*t)))./(c(3)*d(2) - c(2)*d(3) - c(1)*d(4) + d(4)*s + c(4)*(d(1) - t)).^3;   
            else
                out = ((-c(4)*d(3) + c(3)*d(4))*(c(3)*d(2) + c(2)*d(3) - c(1)*d(4) + d(4)*s + c(4)*(-d(1) + t)))./((c(4)*d(1) + c(3)*d(2) - c(2)*d(3) - c(1)*d(4) + d(4)*s - c(4)*t).^2 + 4*(c(4)*d(2) - c(2)*d(4)).*(d(3)*(c(1) - s) + c(3)*(-d(1) + t))).^(3/2);

            end
        end  
        
        function out = d2invT13(obj, s, t)
            
            v = obj.v;
            
            % Compute the change of variables:
            M = [1 -1 -1  1;    % (-1, 1)
                 1  1 -1 -1;    % (1, -1)
                 1  1  1  1;    % (1,  1)
                 1 -1  1 -1];

            % Solve 4x4 linear system:
            c = M \ v(:,1);
            d = M \ v(:,2);

            % Fill in Aaron's formulas for the quad transformations:
            A = c(2)*d(4) - c(4)*d(2);
            if ( A == 0 )
                out = (2*c(4)*(c(3)^2*d(2) - c(3)*(c(2)*d(3) + d(4)*(c(1) - s)) + c(4)*d(3)*(c(1) - s)))./(c(3)*d(2) - c(2)*d(3) - c(1)*d(4) + d(4)*s + c(4)*(d(1) - t)).^3;    
            else
                out = (2*(c(4)*d(3) - c(3)*d(4))*(c(2)*c(3) + c(4)*(-c(1) + s)))./((c(4)*d(1) + c(3)*d(2) - c(2)*d(3) - c(1)*d(4) + d(4)*s - c(4)*t).^2 + 4*(c(4)*d(2) - c(2)*d(4)).*(d(3)*(c(1) - s) + c(3)*(-d(1) + t))).^(3/2);

            end
        end              
        
        function out = d2invT21(obj, s, t)
            
            v = obj.v;
            
            % Compute the change of variables:
            M = [1 -1 -1  1;    % (-1, 1)
                 1  1 -1 -1;    % (1, -1)
                 1  1  1  1;    % (1,  1)
                 1 -1  1 -1];

            % Solve 4x4 linear system:
            c = M \ v(:,1);
            d = M \ v(:,2);

            % Fill in Aaron's formulas for the quad transformations:
            A = c(3)*d(4)-c(4)*d(3);        
            if ( A == 0 )
                out = (2*d(4)*(-c(3)*d(2)^2 + c(4)*d(2)*(d(1) - t) + c(2)*(d(3)*d(2) - d(1)*d(4) + d(4)*t)))./(c(2)*d(3) - c(3)*d(2) - c(1)*d(4) + d(4)*s + c(4)*(d(1) - t)).^3;
            else
                out = (2*(-c(4)*d(2) + c(2)*d(4))*(d(3)*d(2) + d(4)*(-d(1) + t)))./((c(4)*d(1) + c(2)*d(3)-c(3)*d(2) - c(1)*d(4) + d(4)*s - c(4)*t).^2 + 4*(c(4)*d(3) - c(3)*d(4))*(d(2)*(c(1) - s) + c(2)*(-d(1) + t))).^(3/2);
            end
        end 
        
        function out = d2invT22(obj, s, t)
            
            v = obj.v;
            
            % Compute the change of variables:
            M = [1 -1 -1  1;    % (-1, 1)
                 1  1 -1 -1;    % (1, -1)
                 1  1  1  1;    % (1,  1)
                 1 -1  1 -1];

            % Solve 4x4 linear system:
            c = M \ v(:,1);
            d = M \ v(:,2);

            % Fill in Aaron's formulas for the quad transformations:
            A = c(3)*d(4)-c(4)*d(3);        
            if ( A == 0 )
                out = (c(2)*d(4)*(-c(2)*d(3) + c(3)*d(2) + c(1)*d(4) - d(4)*s) + c(4)^2*d(2)*(-d(1) + t) + c(4)*(d(2)*(c(3)*d(2) - c(1)*d(4) + d(4)*s) - c(2)*(d(3)*d(2) - d(1)*d(4) + d(4)*t)))./(c(2)*d(3) - c(3)*d(2) - c(1)*d(4) + d(4)*s + c(4)*(d(1) - t)).^3;
            else
                out = -((-c(4)*d(2) + c(2)*d(4))*(c(2)*d(3) + c(3)*d(2) - c(1)*d(4) + d(4)*s + c(4)*(-d(1) + t)))./((c(4)*d(1) + c(2)*d(3) - c(3)*d(2) - c(1)*d(4) + d(4)*s - c(4)*t).^2 + 4*(c(4)*d(3) - c(3)*d(4)).*(d(2)*(c(1) - s) + c(2)*(-d(1) + t))).^(3/2);
            end
        end 
        
        function out = d2invT23(obj, s, t)
            
            v = obj.v;
            
            % Compute the change of variables:
            M = [1 -1 -1  1;    % (-1, 1)
                 1  1 -1 -1;    % (1, -1)
                 1  1  1  1;    % (1,  1)
                 1 -1  1 -1];

            % Solve 4x4 linear system:
            c = M \ v(:,1);
            d = M \ v(:,2);

            % Fill in Aaron's formulas for the quad transformations:
            A = c(3)*d(4)-c(4)*d(3);        
            if ( A == 0 )
                out = (2*c(4)*(c(2)^2*d(3) - c(2)*(c(3)*d(2) + d(4)*(c(1) - s)) + c(4)*d(2)*(c(1) - s)))./(c(2)*d(3) - c(3)*d(2) - c(1)*d(4) + d(4)*s + c(4)*(d(1) - t)).^3;
            else
                out = (2*(-c(4)*d(2) + c(2)*d(4))*(c(3)*c(2) + c(4)*(-c(1) + s)))./((c(4)*d(1) + c(2)*d(3) - c(3)*d(2) - c(1)*d(4) + d(4)*s - c(4)*t).^2 + 4*(c(4)*d(3) - c(3)*d(4)).*(d(2)*(c(1) - s) + c(2)*(-d(1) + t))).^(3/2);
            end
        end         
        
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STATIC METHODS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods ( Access = public, Static = true )
        
        function v = assertIsQuad(v)
        %ASSERTISQUAD    Check we have valid vertices for a ultraSEMQuad.
        
            if ( ~isnumeric(v) )
                error('Input should be numeric.');
            end
            % Ensure v is of the form [x, y], not [x ; y]:
            if ( size(v,2) ~= 2 ), v = v.'; end
            % Check dimension:
            if ( size(v,2) ~= 2 || size(v, 1) ~= 4 )
                error('Incorrect vertices dimension.')
            end
            % Ensure vertices are oriented in an anticlockwise direction:
            if ( ultraSEMDomain.isClockwise(v) )
                % Switch the second and fourth indices.
                v([2,4],:) = v([4,2],:);
            end
            
        end

    end

end
