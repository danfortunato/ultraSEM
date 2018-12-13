classdef ultraSEMPDO
%ULTRASEMPDO   Partial differential operator (PDO) from the ULTRASEM system.
%   An ULTRASEMPDO object stores the coefficients of a linear PDO of the
%   form
%       Lu = dxx*u_xx + dxy*u_xy + dyy*u_yy + dx*u_x + dy*u_y + b*u.
%
%   The possible syntaxes for constructing an ULTRASEMPDO are
%    * ULTRASEMPDO({dxx, dxy, dyy}, {dx, dy}, b)
%    * ULTRASEMPDO({dxx, dyy}, {dx, dy}, b) - in which case d_xy = 0.
%    * ULTRASEMPDO(dxx, {dx, dy}, b) - in which case d_yy = d_xx and d_xy = 0.
%    * ULTRASEMPDO({dxx, dxy, dyy}, dx, b) - in which case d_y = dx.
%   In all cases, each coefficient may be a scalar or a function handle of
%   the form @(x,y).
%
%   Calls of the form ULTRASEM({{dxx, dxy, dxyy}, {dx, dy}, b}) are also
%   supported.

% Copyright 2018 by Nick Hale and Dan Fortunato.

    %#ok<*PROP>

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties ( Access = public )

        dxx = 0 % Coefficient of (d^2/dx^2).
        dxy = 0 % Coefficient of (d^2/dxdy).
        dyy = 0 % Coefficient of (d^2/dy^2).
        dx  = 0 % Coefficient of (d/dx).
        dy  = 0 % Coefficient of (d/dy).
        b   = 0 % Coefficient of constant term.

        % i.e.,
        %   Lu = dxx*u_xx + dxy*u_xy + dyy*u_yy + dx*u_x + dy*u_y + b*u.

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods

        function obj = ultraSEMPDO(varargin)

            if ( nargin == 0 )
                % Construct empty PDO:
                return
            elseif ( nargin == 3 )
                % Convert to a cell:
                op = varargin;
            elseif ( nargin == 1 && iscell(varargin{1}) )
                op = varargin{1};
            else
                error('Unknown input sequence.')
            end

            % Second-order terms:
            if ( iscell(op{1}) )
                if ( numel(op{1}) == 3 )
                    obj.dxx = op{1}{1};
                    obj.dxy = op{1}{2};
                    obj.dyy = op{1}{3};
                elseif ( numel(op{1}) == 2 )
                    obj.dxx = op{1}{1};
                    obj.dyy = op{1}{2};
                elseif ( numel(op{1}) == 1)
                    obj.dxx = op{1}{1};
                    obj.dyy = op{1}{1};
                end
            elseif ( isnumeric(op{1}) )
                if ( ~isscalar(op{1}) )
                    error('ULTRASEM:ULTRASEMPDO:nonScalar', ...
                        'Numeric input must be a scalar');
                end
                obj.dxx = op{1};
                obj.dyy = op{1};
            elseif ( isa(op{1}, 'function_handle') )
                obj.dxx = op{1};
                obj.dyy = op{1};
            end

            % First-order (advection) terms:
            if ( iscell(op{2}) )
                if ( numel(op{2}) == 2 )
                    obj.dx = op{2}{1};
                    obj.dy = op{2}{2};
                elseif ( numel(op{1}) == 1)
                    % TODO: Do we want this?
                    obj.dx = op{2}{1};
                    obj.dy = op{2}{1};
                end
            elseif ( isnumeric(op{2}) )
               if ( isempty(op{2}) )
                   op{2} = 0;
               elseif ( ~isscalar(op{2}) )
                    error('ULTRASEM:ULTRASEMPDO:nonScalar', ...
                       'Numeric input must be a scalar.');
               end
               % TODO: Do we want this?
               obj.dx = op{2};
               obj.dy = op{2};
            end

            % Constant (reaction) terms:
            obj.b = op{3};

        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% METHODS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods ( Access = public, Static = false )

        function C = pdo2cell(L)
        %PDO2CELL   Convert a PDO to a cell array of coefficients.
        %   C = PDO2CELL(L) extracts the coefficients of the PDO L and
        %   places them in a cell array C = {{dxx, dxy, dyy}, {dx, dy}, b}.
            C = cell(1, 3);
            C{1} = {L.dxx, L.dxy, L.dyy};
            C{2} = {L.dx, L.dy};
            C{3} = L.b;
        end

        [op, isConstant] = feval(op, x, y);

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
