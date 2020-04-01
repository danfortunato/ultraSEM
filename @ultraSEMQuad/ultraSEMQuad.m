classdef ultraSEMQuad < ultraSEMMapping
%ULTRASEMQUAD   Quadrilateral mapping object from ULTRASEM system.

    %#ok<*PROP>

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
        scl % Scale factor for RHS
    end

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

            if ( ~iscell(v) )
                v = {v};
            end
            nobj = numel(v);
            obj = repelem(obj,nobj,1);
            for k = 1:nobj
                v{k} = ultraSEMQuad.assertIsQuad(v{k});
                obj(k).v = v{k};
                obj(k).scl = @(r,s) obj(k).det(r,s).^3;
            end

        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% METHODS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods
        
        function obj = parametrize(obj)
            % Compute the change of variables:
            M = [1 -1 -1  1;  % (-1, 1)
                 1  1 -1 -1;  % (1, -1)
                 1  1  1  1;  % (1,  1)
                 1 -1  1 -1];
            params = M \ obj.v;
            obj.a1 = params(1,1); obj.a2 = params(1,2);
            obj.b1 = params(2,1); obj.b2 = params(2,2);
            obj.c1 = params(3,1); obj.c2 = params(3,2);
            obj.d1 = params(4,1); obj.d2 = params(4,2);
        end

        function c = centroid(Q)
            [xmid, ymid] = centroid(polyshape(Q.v));
            c = [xmid ; ymid]; 
        end
        
        function out = isRect(M)
            out = false;
        end
                
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
            if ( A == 0 )
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
            if ( A == 0 )
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
