<<<<<<< Updated upstream
classdef ultraSEMMapping < handle
=======
classdef ultraSEMMapping < matlab.mixin.Heterogeneous
>>>>>>> Stashed changes
%ULTRASEMMAPPING  Abstract mapping object from ULTRASEM system.

    %#ok<*PROP>

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties ( Access = public )

        T1             % square->domain map 1st coordinate
        T2             % square->domain map 2nd coordinate
        invT1          % domain->square map 1st coordinate
        invT2          % domain->square map 2nd coordinate
        dinvT1         % 1st derivatives of T1^{-1}
        dinvT2         % 1st derivatives of T2^{-1}
        d2invT1        % 2nd derivatives of T1^{-1}
        d2invT2        % 2nd derivatives of T2^{-1}

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false, Abstract = true )

    end

    methods ( Access = public, Static = false, Sealed )

        function [xx, yy] = chebGrid(n, map)
            [xx, yy] = util.chebpts2(n, n, map.domain);
            [xx, yy] = transformGrid(map, xx, yy);
        end

        function n = normals( Q )
        %NORMALS   Outward pointing normal vectors to the edges of a quad.
            v = [ Q.T1([-1 1 1 -1],[-1 -1 1 1]) ;
                  Q.T2([-1 1 1 -1],[-1 -1 1 1]) ].';
            v = [ v ; v(1,:) ];
            dx = diff(v(:,1));
            dy = diff(v(:,2));
            n = [dy.' ; -dx.' ];
            n = n(:,[4 2 1 3]);  % Reorder to left, right, down, up
            n = n * diag(1./sqrt(sum( n.^2 )));  % Normalize
        end

        function plot( T, varargin )
            % Plot the mapped domain and grid:
            
            % Choose a color:
            if ( nargin > 1 && ischar(varargin{1}) && ...
                    ~isempty(regexp( varargin{1}, '[bgrcmykw]', 'match')) )
                col = varargin{1};
                varargin(1) = [];
            elseif ( nargin > 1 && isnumeric(varargin{1}) )
                col = varargin{1};
                varargin(1) = [];
            else
%                 col = 'm';
                col = rand(1, 3);
            end

            plotPts = false;
            n = 21;

            holdState = ishold();

            if ( nargin > 1 && ...
                    ~isempty(regexp( varargin{1}, '[.ox+*sdv^<>ph]', 'match')) )
                plotPts = true;
            end

            for k = 1:numel(T)
                % Plot domain:
                v = vertices(T(k));
                h(k,1) = fill(v(1,[1:end, 1]), v(2,[1:end, 1]), col, ...
                    'FaceAlpha', .25, varargin{:}); hold on %#ok<AGROW>
                % Add text to center of patch:

                if ( plotPts )
                    % Plot the grid:
                    [xx, yy] = chebGrid( n, T(k) );
                    plot(xx, yy, varargin{:})
                end

                % Add text to centre of patch:
                if ( numel(T) < 100 ) % Don't add text on large meshes
                    c = centroid(T(k));
                    text(c(1), c(2), int2str(k), 'HorizontalAlignment', 'center')
                end

            end

            if ( ~holdState )
                hold off
            end

        end
        
        function T = plus(T, c)
        %PLUS   Shift an ultraSEMDomain.
        %   T + C will shift the ultraSEMDomain to the right by real(c) and
        %   upwards by imag(C). C must be a scalar.
        %
        % See also MINUS().

            if ( ~isnumeric(c) )
                error('ULTRASEM:ULTRASEMDOMAIN:plus:unknown', ...
                    'Cannot add an object of type %s to a ultraSEMDomain.', ...
                    class(c));
            elseif ( ~isscalar(c) )
                error('ULTRASEM:ULTRASEMDOMAIN:plus:scalar', ...
                    'C must be a scalar.')
            end

            % Shift the domain:
            for k = 1:numel(T)
                T(k).v(:,1) = T(k).v(:,1) + real(c);
                T(k).v(:,2) = T(k).v(:,2) + imag(c);
            end

        end

        function normal_d = transformNormalD( T, xy, n )
        %TRANSFORMNORMALD   Normal derivative operator for mapped domains.

            % Chebyshev differentiation matrices in (s,t), i.e. mapped
            % [-1,1]^2 space
            I = speye( n );
            D = util.diffmat( n, 1 );
            S = util.convertmat( n, 0, 0 );
            Dcheb = S\D;
            Ds = kron( Dcheb, I );
            Dt = kron( I, Dcheb );

            % Evaluation matrices on the four sides
            lbc = kron( (-1).^(0:n-1), I );
            rbc = kron( ones(1,n), I );
            dbc = kron( I, (-1).^(0:n-1) );
            ubc = kron( I, ones(1,n) );

            % (s,t)-differentiation matrices on the four sides
            Ds_l = lbc * Ds; Dt_l = lbc * Dt;
            Ds_r = rbc * Ds; Dt_r = rbc * Dt;
            Ds_d = dbc * Ds; Dt_d = dbc * Dt;
            Ds_u = ubc * Ds; Dt_u = ubc * Dt;

            % Derivatives of (s,t) with respect to (x,y)
            dsdx = @(x,y) T.dinvT1{1}(x,y);
            dsdy = @(x,y) T.dinvT1{2}(x,y);
            dtdx = @(x,y) T.dinvT2{1}(x,y);
            dtdy = @(x,y) T.dinvT2{2}(x,y);

            % Boundary points on the four sides
            xxl = xy{1}(:,1); yyl = xy{1}(:,2);
            xxr = xy{2}(:,1); yyr = xy{2}(:,2);
            xxd = xy{3}(:,1); yyd = xy{3}(:,2);
            xxu = xy{4}(:,1); yyu = xy{4}(:,2);

            % Get the 1D Chebyshev coefficients of d(s,t)/d(x,y) on the
            % four sides through evaluation at the boundary nodes, and
            % construct multiplication matrices
            coeffMult = @(vals) util.multmat(n, util.vals2coeffs(vals), 0);
            dsdx_l = coeffMult( dsdx(xxl,yyl) ); dsdy_l = coeffMult( dsdy(xxl,yyl) );
            dtdx_l = coeffMult( dtdx(xxl,yyl) ); dtdy_l = coeffMult( dtdy(xxl,yyl) );
            dsdx_r = coeffMult( dsdx(xxr,yyr) ); dsdy_r = coeffMult( dsdy(xxr,yyr) );
            dtdx_r = coeffMult( dtdx(xxr,yyr) ); dtdy_r = coeffMult( dtdy(xxr,yyr) );
            dsdx_d = coeffMult( dsdx(xxd,yyd) ); dsdy_d = coeffMult( dsdy(xxd,yyd) );
            dtdx_d = coeffMult( dtdx(xxd,yyd) ); dtdy_d = coeffMult( dtdy(xxd,yyd) );
            dsdx_u = coeffMult( dsdx(xxu,yyu) ); dsdy_u = coeffMult( dsdy(xxu,yyu) );
            dtdx_u = coeffMult( dtdx(xxu,yyu) ); dtdy_u = coeffMult( dtdy(xxu,yyu) );

            % (x,y)-differentiation matrices on the four sides, computed
            % via the chain rule
            Dx_l = dsdx_l * Ds_l + dtdx_l * Dt_l;
            Dy_l = dsdy_l * Ds_l + dtdy_l * Dt_l;
            Dx_r = dsdx_r * Ds_r + dtdx_r * Dt_r;
            Dy_r = dsdy_r * Ds_r + dtdy_r * Dt_r;
            Dx_d = dsdx_d * Ds_d + dtdx_d * Dt_d;
            Dy_d = dsdy_d * Ds_d + dtdy_d * Dt_d;
            Dx_u = dsdx_u * Ds_u + dtdx_u * Dt_u;
            Dy_u = dsdy_u * Ds_u + dtdy_u * Dt_u;

            % Convert to normal derivatives on each side
            vn = normals(T);
            normal_d = [ vn(1,1)*Dx_l + vn(2,1)*Dy_l ;  % "Left"
                         vn(1,2)*Dx_r + vn(2,2)*Dy_r ;  % "Right"
                         vn(1,3)*Dx_d + vn(2,3)*Dy_d ;  % "Down"
                         vn(1,4)*Dx_u + vn(2,4)*Dy_u ]; % "Up"

        end

        function [X, Y] = transformGrid( T, x, y )
            X = T.T1(x, y);
            Y = T.T2(x, y);
        end
        
        function v = vertices(M)
%             xRef = M.domain([1 2 2 1]);
%             yRef = M.domain([3 3 4 4]);
            xRef = [-1 1 1 -1]; yRef = [-1 -1 1 1];
            x = M.T1(xRef, yRef);
            y = M.T2(xRef, yRef);
            v = [x ; y];
        end

    end

    methods (Static, Sealed, Access = protected)
        function cobj = convertObject(~, v)
            n = size(v, 1);
            cobj(n,1) = rectangle();
            for k = 1:n
                cobj(k,1) = rectangle(v(k,:));
            end
        end
    end

end
