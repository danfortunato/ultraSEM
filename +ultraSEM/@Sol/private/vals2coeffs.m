function C = vals2coeffs(V)
%VALS2COEFFS   Convert a cell array of values to 2D Chebyshev coefficients.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

    C = cell(size(V));
    for k = 1:length(V)
        C{k} = util.vals2coeffs(util.vals2coeffs(V{k}).').';
    end
end
