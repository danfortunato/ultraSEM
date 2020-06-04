function coeffs = bc2coeffs(S, bc)
%BC2COEFFS   Convert boundary conditions to Chebyshev coefficients.

% TODO - Document.

assert(isBuilt(S), 'ultraSEM object must be built before creating BCs.');

% Get boundary points:
edges = S.patches{1}.edges; % TODO: Breaks encapsulation

% Initialize boundary data. We should have n coefficients for
% each boundary.
coeffs = cell(size(edges));
if ( ~isnumeric(bc) )
    % Evaluate the BC if given a function handle:
    for k = 1:size(coeffs, 1)
        % Create grid on this edge:
        a = edges(k,1)+edges(k,2)*1i;
        b = edges(k,3)+edges(k,4)*1i;
        t = chebpts(edges(k,5));
        x = real(b-a)/2*t+real(b+a)/2;
        y = imag(b-a)/2*t+imag(b+a)/2;
        % Evaluate at grid:
        vals = feval(bc, x, y);
        % Convert from values to coeffs:
        coeffs{k} = util.vals2coeffs(vals);
    end
elseif ( isscalar(bc) )
    % Convert a scalar to coeffs:
    for k = 1:size(coeffs, 1)
        n = edges(k,5);
        coeffs{k} = [bc ; zeros(n-1, 1)];
    end
else
    error('ULTRASEM:ULTRASEM:bc2coeffs:badBC', ...
        'Cannot evaluate boundary data.');
end

% CAT() is 10x faster than CELL2MAT().
coeffs = cat(1, coeffs{:});

end
