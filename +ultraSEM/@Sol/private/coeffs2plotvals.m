function V = coeffs2plotvals(C)
%COEFFS2PLOTVALS   Convert 2D Chebyshev coefficients to values at plot points.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

persistent Eval
nplotpts = ultraSEM.Sol.nplotpts;
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
