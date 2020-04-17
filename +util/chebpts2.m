function [xx, yy] = chebpts2(nx, ny, dom)
%CHEBPTS2   Chebyshev tensor product grid.
%   [XX, YY] = CHEBPTS2(N) constructs an N x N grid of tensor-product
%   Chebyshev points on [-1,1]^2.
%
%   [XX, YY] = CHEBPTS2(NX, NY) constructs an NX x NY grid of
%   tensor-product Chebyshev points on [-1,1]^2.
%
%   [XX, YY] = CHEBPTS2(NX, NY, DOM) constructs an NX x NY grid of
%   tensor-product Chebyshev points on the rectangle [a,b] x [c,d], where
%   DOM = [a b c d]. If DOM is an ULTRASEM.MAPPING, then the points are
%   mapped from [-1,1]^2 according to DOM.
%
% See also CHEBPTS.

if ( nargin == 1 )
    ny = nx;
end

map = [];
if ( nargin < 3 )
    dom = [-1 1 -1 1];
elseif ( isa(dom, 'ultraSEM.Mapping') )
    map = dom;
    dom = [-1 1 -1 1];
end

x = util.chebpts(nx, dom(1:2));
y = util.chebpts(ny, dom(3:4));

% Tensor product
[xx, yy] = meshgrid(x, y);

if ( ~isempty(map) )
    [xx, yy] = transformGrid(map, xx, yy);
end

end
