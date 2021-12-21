function I = integral2(f, varargin)
%INTEGRAL2   Double integral of an ULTRASEM.SOL.
%   I = INTEGRAL2(F) returns the double integral of the ULTRASEM.SOL F over
%   its domain.
%
%   I = INTEGRAL2(F, 'all') returns an array of double integrals over each
%   patch of F.
%
%   See also SUM2.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

% Parse arguments.
reduce = true;
if ( nargin == 2 )
    if ( strcmp(varargin{1}, 'all') )
        reduce = false;
    end
end

% Empty check:
if ( isempty(f) )
    I = 0;
    return
end

% If a patch uses an N x N discretization, then quadrature is performed on
% that patch using N points.
I = zeros(length(f),1);
for k = 1:length(f)
    [ny,nx] = size(f.coeffs{k});
    qx = nx; qy = ny;
    U = zeros(qy,qx);
    U(1:ny,1:nx) = f.coeffs{k}; U = U(1:qy,1:qx);
    V = util.coeffs2vals( util.coeffs2vals(U).' ).';
    wx = util.quadwts(qx); wx = wx(:);
    wy = util.quadwts(qy); wy = wy(:);
    [rr, ss] = util.chebpts2(qx, qy);
    jac = f.domain(k).det(rr, ss);
    I(k) = sum(sum(V .* wy .* wx.' .* jac));
end

% Combine norms on each element.
if ( reduce )
    I = sum(I);
end

end
