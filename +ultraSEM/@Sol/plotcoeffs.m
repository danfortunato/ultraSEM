function plotcoeffs(sol)
%PLOTCOEFFS   Display coefficients graphically.
%   PLOTCOEFFS(SOL) plots the absolute values of the coefficients
%   underlying the representation of the ULTRASEM.SOL SOL on a semilog
%   scale.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

holdState = ishold();

% Loop over the patches:
[x,y] = getGrid(sol);
for k = 1:length(sol)
    stem3(x{k}, y{k}, abs(sol.coeffs{k}), 'ok', 'MarkerFaceColor', 'k');
    hold on
end
set(gca, 'ZScale', 'log')

if ( ~holdState )
    hold off;
end

end
