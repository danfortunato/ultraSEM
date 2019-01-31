classdef ultraSEMSol
%ULTRASEMSOL  Solution object from ULTRASEM system.
%   An ULTRASEMSOL is the result of solving a PDE using ultraSEM. It allows
%   elementary exploration of the solution, such as various methods of
%   plotting and evaluation within the domain.

% Copyright 2018 by Nick Hale and Dan Fortunato.

    %#ok<*PROP>

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties ( Access = public )

        domain % Domain of the patch.
        x      % Interpolation points on each patch in x.
        y      % Interpolation points on each patch in y.
        u      % Coefficients of solution.

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods

        function obj = ultraSEMSol(u, d)
            %ULTRASEMSOL   Class constructor for the @ultraSEMSol class.

            % Construct empty patch:
            if ( nargin == 0 )
                return
            end

            obj.domain = d;
            obj.u = u;

            for k = 1:size(d, 1)
                n = size(u{k}, 1);
                [obj.x{k,1}, obj.y{k,1}] = chebpts2(n, n, d(k,:));
            end

        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% METHODS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods ( Access = public, Static = false )

        function chebpolyplot(sol)

            holdState = ishold();

            % Loop over the patches:
            for k = 1:length(sol)
                stem3(sol.x{k}, sol.y{k}, abs(sol.u{k}), 'ok', 'MarkerFaceColor', 'k');
                hold on
            end
            set(gca, 'ZScale', 'log')

            if ( ~holdState )
                hold off;
            end
        end

        function chebpolyfill(sol)

            holdState = ishold();
            n = size(sol.u{1}, 1);

            % Loop over the patches:
            data = {};
            for k = 1:length(sol)
                u = abs(sol.u{k});
                tail = [u(end-floor(.1*n):end,:) ; u(:,end-floor(.1*n):end)'];
                err(k) = norm(tail(:),inf)./max(u(:));
                err(k) = log10(err(k));
                d = sol.domain(k,:);
                data = [data, d([1 1 2 2 1]), d([3 4 4 3 3]), err(k)];
            end

            fill(data{:});

            if ( ~holdState )
                hold off;
            end
        end

        function contour(sol, varargin)

            holdState = ishold();
            N = 10;
            levels = [];

            % See if an N was given:
            if ( nargin > 1 && isnumeric(varargin{1}) )
                v1 = varargin{1};
                if ( isscalar(v1) )
                    N = v1;
                else
                    levels = v1;
                end
                varargin(1) = [];
            end

            % Determine some levels:
            if ( isempty(levels) )
                minu = minEst(sol);
                maxu = maxEst(sol);
                levels = linspace(minu, maxu, N);
            end

            % Loop over the patches:
            vals = coeffs2vals(sol.u);
            for k = 1:length(sol)
                u = vals{k};
                if ( ~isreal(u) )
                    u = abs(u);
                end
                contour(sol.x{k}, sol.y{k}, u, levels, varargin{:});
                hold on
            end
            view(0, 90)

            % Fix limits (we don't want to do axis equal because of z)
            xl = xlim; yl = ylim;
            if ( diff(xl) < diff(yl) )
                xlim(xl+(diff(yl)-diff(xl))*[-.5 .5])
            else
                ylim(yl+(diff(xl)-diff(yl))*[-.5 .5])
            end

            if ( ~holdState )
                hold off
            end

        end

        function u = feval(sol, x, y)

            if ( ~isnumeric(sol.domain) )
                u = fevalMapped(sol, x, y);
                return
            end

            minx = min(sol.domain(:,1));
            miny = min(sol.domain(:,3));

            sx = size(x);
            u = 0*x(:);
            x = x(:);
            y = y(:);

            for k = 1:length(sol)

                dom = sol.domain(k,:);
                idx = (x > dom(1) | (dom(1) == minx & x == minx)) & x <= dom(2) & ...
                    (y > dom(3) | (dom(3) == miny & y == miny)) & y <= dom(4);

                if ( any(idx) )
                    % Map x and y to [-1,1]^2:
                    sclx = 2/diff(dom(1:2));
                    scly = 2/diff(dom(3:4));
                    xm = sclx*(x(idx)-dom(1))-1;
                    ym = scly*(y(idx)-dom(3))-1;
                    u(idx) = clenshaw(sol.u{k}, xm, ym);
                end

            end

            u = reshape(u, sx);

        end

        function u = fevalMapped(sol, x, y)

            sx = size(x);
            u = 0*x(:);
            x0 = x(:);
            y0 = y(:);

            for k = 1:size(sol.domain, 1 )

                map = sol.domain(k,:);
                x = map.invT1(x0, y0);
                y = map.invT2(x0, y0);
                dom = map.domain;

                idx = ( isreal(x) & isreal(y) & x >= dom(1) & x <= dom(2) & ...
                       y >= dom(3) & y <= dom(4));

                if ( any(idx) )
                    u(idx) = clenshaw(sol.u{k}, x(idx), y(idx));
                end

            end

            u = reshape(u, sx);

        end

        function out = length(sol)
            out = length(sol.u);
        end

        function mesh(sol, varargin)

            holdState = ishold();

            vals = coeffs2vals(sol.u);
            for k = 1:length(sol)
                u = vals{k};
                if ( ~isreal(u) )
                    u = abs(u);
                end
                mesh(sol.x{k}, sol.y{k}, u, varargin{:});
                hold on
            end
            view(0, 90)

            % Fix limits (we don't want to do axis equal because of z)
            xl = xlim; yl = ylim;
            if ( diff(xl) < diff(yl) )
                xlim(xl+(diff(yl)-diff(xl))*[-.5 .5])
            else
                ylim(yl+(diff(xl)-diff(yl))*[-.5 .5])
            end

            if ( ~holdState )
                hold off
            end

        end

        function m = maxEst(sol)
            m = max(cellfun(@(u) max(u(:)), coeffs2vals(sol.u)));
        end

        function m = minEst(sol)
            m = min(cellfun(@(u) min(u(:)), coeffs2vals(sol.u)));
        end

        function h = minus( f, g )
            h = plus( f, (-g) );
        end

        function val = normest( sol )
           % Estimate maximum absolute value of F over all patches.
            val = max(cellfun(@(u) max(abs(u(:))), coeffs2vals(sol.u)));
        end

        function varargout = plot(varargin)

            [varargout{1:nargout}] = surf(varargin{:});

        end

        function h = plus( f, g )
            % PLUS for ultraSEMSol

            if ( isnumeric( f ) )
                h = plus(g, f);
                return
            elseif ( isnumeric( g ) )
                h = f;
                for k = 1:size(h.u,1)
                    h.u{k}(1,1) = h.u{k}(1,1) + g;
                end
            elseif ( isa(f, 'ultraSEMSol') && isa(g, 'ultraSEMSol') )
                h = f;
                % TODO: Assume on the same grid for now.
                h.u = cellfun(@plus, f.u , g.u, 'UniformOutput', false);
            end
        end

        function scatter(sol, varargin)

            numPatches = size(sol.domain, 1);
            holdState = ishold();

            vals = coeffs2vals(sol.u);
            for k = 1:numPatches
                scatter(sol.x{k}, sol.y{k}, vals{k}, 'edgealpha', 0, varargin{:});
                hold on
            end
            view(0, 90)

            % Fix limits (we don't want to do axis equal because of z)
            xl = xlim; yl = ylim;
            if ( diff(xl) < diff(yl) )
                xlim(xl+(diff(yl)-diff(xl))*[-.5 .5])
            else
                ylim(yl+(diff(xl)-diff(yl))*[-.5 .5])
            end

            if ( ~holdState )
                hold off
            end

        end

        function varargout = surf(sol, varargin)

            holdState = ishold();

            vals = coeffs2vals(sol.u);
            for k = 1:length(sol)
                u = vals{k};
                if ( ~isreal(u) )
                    u = abs(u);
                end
                h(k) = surf(sol.x{k}, sol.y{k}, u, 'edgealpha', 1, varargin{:});
                hold on
            end
            shading interp
            view(0, 90)

            % Fix limits (we don't want to do axis equal because of z)
            xl = xlim; yl = ylim;
            if ( diff(xl) < diff(yl) )
                xlim(xl+(diff(yl)-diff(xl))*[-.5 .5])
            else
                ylim(yl+(diff(xl)-diff(yl))*[-.5 .5])
            end

            if ( ~holdState )
                hold off
            end

            if ( nargout > 0 )
                varargout = {h};
            end

        end

        function varargout = pcolor(sol, varargin)

            holdState = ishold();

            vals = coeffs2vals(sol.u);
            for k = 1:length(sol)
                u = vals{k};
                if ( ~isreal(u) )
                    u = abs(u);
                end
                h(k) = surface(sol.x{k}, sol.y{k}, 0*u, u, 'edgealpha', 1, varargin{:});
                hold on
            end
            shading interp
            view(0, 90)

            % Fix limits (we don't want to do axis equal because of z)
            xl = xlim; yl = ylim;
            if ( diff(xl) < diff(yl) )
                xlim(xl+(diff(yl)-diff(xl))*[-.5 .5])
            else
                ylim(yl+(diff(xl)-diff(yl))*[-.5 .5])
            end

            if ( ~holdState )
                hold off
            end

            if ( nargout > 0 )
                varargout = {h};
            end

        end

        function f = uplus( f )
        end

        function f = uminus( f )
            f.u = cellfun(@uminus, f.u, 'UniformOutput', false);
        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STATIC METHODS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods ( Access = public, Static = true )

    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% METHODS IMPLEMENTED IN THIS FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function V = coeffs2vals(C)
%COEFFS2VALS   Convert a cell array of 2D Chebyshev coefficients to values.
    V = cell(size(C));
    for k = 1:length(C)
        V{k} = chebtech2.coeffs2vals(chebtech2.coeffs2vals(C{k}).').';
    end
end

function v = clenshaw(C, x, y)
%CLENSHAW   Evaluate a 2D Chebyshev expansion at the given points.
    v = 0*x;
    Cy = chebtech2.clenshaw(y, C).';
    for k = 1:numel(x)
        v(k) = chebtech2.clenshaw(x(k), Cy(:,k));
    end
end
