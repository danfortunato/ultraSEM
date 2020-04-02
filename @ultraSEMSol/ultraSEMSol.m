classdef ultraSEMSol
%ULTRASEMSOL  Solution object from ULTRASEM system.
%   An ULTRASEMSOL is the result of solving a PDE using ultraSEM. It allows
%   elementary exploration of the solution, such as various methods of
%   plotting and evaluation within the domain.

% Copyright 2018 by Nick Hale and Dan Fortunato.

    %#ok<*PROP>
    %#ok<*PROPLC>

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties ( Access = public )

        domain % Domain of the patch.
        u      % Coefficients of solution.

    end

    properties (Constant)

        nplotpts = 200; % Number of plot points.

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

        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% METHODS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods ( Access = public, Static = false )

        function chebpolyplot(sol)

            holdState = ishold();

            % Loop over the patches:
            [x,y] = getGrid(sol);
            for k = 1:length(sol)
                stem3(x{k}, y{k}, abs(sol.u{k}), 'ok', 'MarkerFaceColor', 'k');
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
            [x, y] = plotpts(sol);
            for k = 1:length(sol)
                u = coeffs2plotvals(sol.u{k});
                if ( ~isreal(u) )
                    u = abs(u);
                end
                contour(x{k}, y{k}, u, levels, varargin{:});
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
                    u(idx) = util.clenshaw2d(sol.u{k}, xm, ym);
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
                % Convert to single to avoid mapping issues
                x = map.r(x0, y0);
                y = map.s(x0, y0);
                dom = [-1 1 -1 1];

                xs = single(x);
                ys = single(y);
                
                idx = ( isreal(xs) & isreal(ys) & xs >= dom(1) & xs <= dom(2) & ...
                       ys >= dom(3) & ys <= dom(4) );

                if ( any(idx) )
                    u(idx) = util.clenshaw2d(sol.u{k}, x(idx), y(idx));
                end

            end

            u = reshape(u, sx);

        end

        function [x, y] = getGrid(sol, kk)

            d = sol.domain;
            u = sol.u;

            x = cell(size(u));
            y = cell(size(u));

            if ( nargin == 1 )
                kk = 1:size(d, 1);
            end

            for k = kk
                nk = size(u{k}, 1);
                if ( isnumeric(d(k,:)) )
                    [x{k,1}, y{k,1}] = util.chebpts2(nk, nk, d(k,:));
                else
                    [xk, yk] = chebpts2(nk);
                    [x{k,1}, y{k,1}] = transformGrid(d(k,:), xk, yk);
                end
            end

            if ( nargin > 1 && numel(kk) == 1 )
                x = x{1};
                y = y{1};
            end

        end

        function [x, y] = plotpts(sol, kk)

            d = sol.domain;
            u = sol.u;

            x = cell(size(u));
            y = cell(size(u));

            if ( nargin < 2 )
                kk = 1:size(d, 1);
            end

            for k = kk
                if ( isnumeric(d(k,:)) )
                    [x{k,1}, y{k,1}] = meshgrid(linspace(d(k,1), d(k,2), sol.nplotpts), ...
                                                linspace(d(k,3), d(k,4), sol.nplotpts));
                else
                    [xk, yk] = meshgrid(linspace(-1, 1, sol.nplotpts));
                    [x{k,1}, y{k,1}] = transformGrid(d(k,:), xk, yk);
                end
            end

            if ( nargin > 1 && numel(kk) == 1 )
                x = x{1};
                y = y{1};
            end

        end

        function out = length(sol)
            out = length(sol.u);
        end

        function mesh(sol, varargin)

            holdState = ishold();

            vals = coeffs2vals(sol.u);
            [x,y] = getGrid(sol);
            for k = 1:length(sol)
                u = vals{k};
                if ( ~isreal(u) )
                    u = abs(u);
                end
                mesh(x{k}, y{k}, u, varargin{:});
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
            [x,y] = getGrid(sol);
            for k = 1:numPatches
                scatter(x{k}, y{k}, vals{k}, 'EdgeAlpha', 0, varargin{:});
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
        
        function varargout = subsref(u, index)
            
            idx = index(1).subs;
            switch index(1).type

            %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% FEVAL / COMPOSE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                case '()'

                    % Where to evaluate:
                    x = idx{1}; 
                    y = idx{2}; 

                    out = feval(u, x, y);

            %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GET %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                case '.'

                    % Call GET() for .PROP access.
                    out = u.(idx);

                otherwise

                    error('CHEBFUN:CHEBFUN:subsref:unexpectedType',...
                        ['??? Unexpected index.type of ', index(1).type]);
            end

            % Convert to a cell:
            varargout = {out};

        end

        function varargout = surf(sol, varargin)

            holdState = ishold();

            [x, y] = plotpts(sol);
            if ( ~iscell(x) )
                x = {x};
                y = {y};
            end
            
            for k = 1:length(sol)
                u = coeffs2plotvals(sol.u{k});
                if ( ~isreal(u) )
                    u = abs(u);
                end
                h(k) = surf(x{k}, y{k}, u, 'EdgeAlpha', 1, varargin{:});
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

            [x,y] = plotpts(sol);

            for k = 1:length(sol)
                u = coeffs2plotvals(sol.u{k});
                if ( ~isreal(u) )
                    u = abs(u);
                end
                h(k) = surface(x{k}, y{k}, 0*u, u, 'EdgeAlpha', 1, varargin{:});
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
        V{k} = util.coeffs2vals(util.coeffs2vals(C{k}).').';
    end
end

function V = coeffs2plotvals(C)

    persistent Eval
    nplotpts = ultraSEMSol.nplotpts;
    pmax = 100; % Maximum p to precompute

    if ( isempty(Eval) )
        Eval = zeros(nplotpts, pmax);
        x = linspace(-1, 1, nplotpts).';
        c = zeros(pmax, 1);
        for k = 1:nplotpts
            c(k) = 1;
            Eval(:,k) = util.clenshaw(c, x);
            c(k) = 0;
        end
    end

    [px,py] = size(C);
    V = Eval(:,1:py) * C * Eval(:,1:px)';

end
