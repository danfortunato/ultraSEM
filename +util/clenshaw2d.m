function v = clenshaw2d(C, x, y)
%CLENSHAW2D   Evaluate a 2D Chebyshev expansion at the given points.
%
%   See also CLENSHAW.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

v = 0*x;
Cy = util.clenshaw(C, y).';
for k = 1:numel(x)
    v(k) = util.clenshaw(Cy(:,k), x(k));
end

end
