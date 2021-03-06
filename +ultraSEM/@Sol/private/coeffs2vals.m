function V = coeffs2vals(C)
%COEFFS2VALS   Convert a cell array of 2D Chebyshev coefficients to values.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

    V = cell(size(C));
    for k = 1:length(C)
        V{k} = util.coeffs2vals(util.coeffs2vals(C{k}).').';
    end
end
