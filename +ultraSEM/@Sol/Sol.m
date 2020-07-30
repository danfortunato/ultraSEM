classdef Sol
%ULTRASEM.SOL   Solution object from the ULTRASEM system.
%   An ULTRASEM.SOL is the result of solving a PDE using ULTRASEM. It
%   allows elementary exploration of the solution, such as various methods 
%   of plotting and evaluation within the domain.
%
%   ULTRASEM.SOL(DOM) constructs an ULTRASEM.SOL representing the zero
%   function over the ULTRASEM.DOMAIN DOM. By default, the discretization
%   size is ULTRASEM.PREF().DISCSIZE. ULTRASEM.SOL(N, DOM) uses an N x N
%   discretization on each patch of DOM. If N is a vector of length
%   LENGTH(DOM), then an N(K) x N(K) discretization is used on patch K of
%   DOM.
%
%   ULTRASEM.SOL(FUNC, DOM) constructs an ULTRASEM.SOL representing the
%   bivariate function handle FUNC over the ULTRASEM.DOMAIN DOM. By
%   default, the discretization size is ULTRASEM.PREF().DISCSIZE.
%   ULTRASEM.SOL(FUNC, N, DOM) uses an N x N discretization on each patch
%   of DOM. If N is a vector of length LENGTH(DOM), then an N(K) x N(K)
%   discretization is used on patch K of DOM.
%
%   ULTRASEM.SOL(COEFFS, DOM) constructs an ULTRASEM.SOL from a cell array
%   of bivariate Chebyshev coefficients. The matrix COEFFS{K} represents
%   the bivariate Chebyshev coefficients of a function over patch K of DOM.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties ( Access = public )

        domain % Domain of the patch.
        coeffs % Coefficients of solution.

    end

    properties ( Constant )

        nplotpts = 200; % Number of plot points.

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods

        function obj = Sol(varargin)

            pref = ultraSEM.Pref;
            n = pref.discSize;
            func = @(x,y) 0*x;
            coeffs = {};

            isValidDom = @(dom) isa(dom, 'ultraSEM.Domain') || ...
                                isa(dom, 'ultraSEM.Mapping');
            isValidDisc = @(n, dom) isnumeric(n) && ...
                ( isscalar(n) || (isvector(n) && length(n) == length(dom)) );

            if ( nargin == 0 )
                % Construct empty ULTRASEM.SOL:
                return
            elseif ( nargin == 1 && isValidDom(varargin{1}) )
                % Call is: SOL(DOM)
                dom = varargin{1};
            elseif ( nargin == 2 && isValidDom(varargin{2}) )
                dom = varargin{2};
                if ( isa(varargin{1}, 'function_handle') )
                    % Call is: SOL(FUNC, DOM)
                    func = varargin{1};
                elseif ( isValidDisc(varargin{1}, dom) )
                    % Call is: SOL(N, DOM)
                    n = varargin{1};
                elseif ( iscell(varargin{1}) && length(varargin{1}) == length(dom) )
                    % Call is: SOL(COEFFS, DOM)
                    coeffs = varargin{1};
                else
                    error('ULTRASEM:SOL:sol:invalid', ...
                        'Invalid call to ultraSEM.Sol constructor.');
                end
            elseif ( nargin == 3 && isa(varargin{1}, 'function_handle') && ...
                     isValidDom(varargin{3}) && isValidDisc(varargin{2}, varargin{3}) )
                % Call is: SOL(FUNC, N, DOM)
                func = varargin{1};
                n = varargin{2};
                dom = varargin{3};
            else
                error('ULTRASEM:SOL:sol:invalid', ...
                    'Invalid call to ultraSEM.Sol constructor.');
            end

            if ( isa(dom, 'ultraSEM.Domain') )
                dom = dom.domain;
            end
            obj.domain = dom;

            if ( isscalar(n) )
                n = repmat(n, length(dom), 1);
            end

            if ( isempty(coeffs) )
                coeffs = cell(length(dom), 1);
                for k = 1:length(dom)
                    if ( isnumeric(dom(k,:)) )
                        [xx, yy] = util.chebpts2(n(k), n(k), dom(k,:));
                    else
                        [xx, yy] = util.chebpts2(n(k));
                        [xx, yy] = transformGrid(dom(k,:), xx, yy);
                    end
                    vals = feval(func, xx, yy);
                    coeffs{k} = util.vals2coeffs( util.vals2coeffs( vals ).' ).';
                end
            end
            obj.coeffs = coeffs;

        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods

        function n = numArgumentsFromSubscript(obj,s,indexingContext)
        %NUMARGUMENTSFROMSUBSCRIPT   Number of arguments for customized indexing methods.
        %   Overloading NUMEL() gives the wrong NARGOUT for SUBSREF().
        %   Defining this function fixes it.
        %
        % See also NUMEL, NARGOUT, SUBSREF.

            n = 1;

        end

    end

end
