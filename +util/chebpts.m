function x = chebpts(n, dom)
%CHEBPTS   Chebyshev points.
%   X = CHEBPTS(N) returns N Chebyshev points on [-1,1].
%
%   X = CHEBPTS(N, DOM) returns N Chebyshev points on the interval DOM.
%
%   See also CHEBPTS2.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

% Default domain:
if ( nargin == 1 )
    dom = [-1 1];
end

% Deal with trivial cases:
if ( n == 0 )
    x = [];
elseif ( n == 1 )
    x = 0;
else
    % Enforce symmetry:
    m = n-1;
    x = sin(pi*(-m:2:m)/(2*m)).';
end

if ( ~all(dom == [-1 1]) )
    % Scale the nodes
    x = dom(1)*(1-x)/2 + dom(2)*(x+1)/2;
end

end
