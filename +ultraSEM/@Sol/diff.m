function f = diff(f, n, dim)
%DIFF   Differentiate an ULTRASEM.SOL.
%   DIFF(F, N, DIM) is the N-th derivative of the ULTRASEM.SOL F in the
%   dimension DIM.
%      DIM = 1 (default) is the derivative in the y direction.
%      DIM = 2 is the derivative in the x direction.
%
%   DIFF(F, [NX NY]) takes the NX-th derivative of F in the x direction
%   and the NY-th derivative of F in the y direction.
%
%   See also DIFFX, DIFFY.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

% Empty check:
if ( isempty(f) )
    return
end

% Default to first derivative:
if ( nargin < 2 || isempty(n) )
    n = 1;
elseif ( n == 0 )
    return
end

% Default to partial derivative in y:
if ( nargin < 3 )
    dim = 1;
elseif ( numel(dim) ~= 1 )
    error('ULTRASEM:SOL:diff:dim', 'Dimension should be either 1 or 2.');
end

% Make sure N is in [NX NY] format
if ( isscalar(n) )
    if ( dim == 1 )
        n = [0 n];
    elseif ( dim == 2 )
        n = [n 0];
    else
        error('ULTRASEM:SOL:diff:dim', 'Dimension should be either 1 or 2.');
    end
end
nx = n(1);
ny = n(2);

for k = 1:length(f)
    coeffs = f.coeffs{k};
    dom = f.domain(k);
    for m = 1:ny, coeffs = mappedDiff(coeffs, dom, 1); end
    for m = 1:nx, coeffs = mappedDiff(coeffs, dom, 2); end
    f.coeffs{k} = coeffs;
end

end

function u = mappedDiff(u, dom, dim)
%DIFFONCE   Mapped derivative of a matrix of Chebyshev coefficients.
%   MAPPEDDIFF(U, DOM, DIM) returns the matrix of Chebyshev coefficients
%   representing the mapped derivative of U in the dimension DIM. The
%   mapping is defined by the coordinate transformations in DOM.

    % Compute derivatives on [-1,1]^2
    [ny, nx] = size(u);
    dfdr = [ cdiff(u.').', zeros(ny, 1) ];
    dfds = [ cdiff(u);     zeros(1, nx) ];

    % Convert to values to multiply by Jacobian factors
    dfdr = util.coeffs2vals( util.coeffs2vals(dfdr).' ).';
    dfds = util.coeffs2vals( util.coeffs2vals(dfds).' ).';

    % Get Jacobian factors for the specified dimension
    [rr, ss] = util.chebpts2(nx, ny);
    if ( dim == 1 )
        dr = @(r,s) dom.drdy(r,s) ./ dom.det(r,s);
        ds = @(r,s) dom.dsdy(r,s) ./ dom.det(r,s);
    elseif ( dim == 2 )
        dr = @(r,s) dom.drdx(r,s) ./ dom.det(r,s);
        ds = @(r,s) dom.dsdx(r,s) ./ dom.det(r,s);
    else
        error('ULTRASEM:SOL:mappedDiff:dim',...
            'Dimension should be either 1 or 2.');
    end

    % Compute df/dx = (dr/dx) df/dr + (ds/dx) df/ds
    vals = dr(rr,ss) .* dfdr + ds(rr,ss) .* dfds;

    % Convert back to coefficients
    u = util.vals2coeffs( util.vals2coeffs(vals).' ).';

end

function dC = cdiff(C)
%CDIFF   Recurrence relation for coefficients of derivative.
%   CDIFF(C) returns the matrix of Chebyshev coefficients whose columns are
%   the derivatives of the columns of C.

    [n, m] = size(C);
    dC = zeros(n-1, m);                        % Initialize vector {c_r}
    w = repmat(2*(1:n-1)', 1, m);
    v = w.*C(2:end,:);                         % Temporal vector
    dC(n-1:-2:1,:) = cumsum(v(n-1:-2:1,:), 1); % Compute c_{n-2}, c_{n-4}, ...
    dC(n-2:-2:1,:) = cumsum(v(n-2:-2:1,:), 1); % Compute c_{n-3}, c_{n-5}, ...
    dC(1,:) = .5*dC(1,:);                      % Adjust the value for c_0

end
