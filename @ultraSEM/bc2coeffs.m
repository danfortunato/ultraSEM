function coeffs = bc2coeffs(S, bc)
%BC2COEFFS   Convert boundary conditions to Chebyshev coefficients.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

% TODO - Document.

assert(isBuilt(S), 'ultraSEM object must be built before creating BCs.');

% Get boundary points:
edges = S.patches{1}.edges; % TODO: Breaks encapsulation

% Initialize boundary data:
coeffs = cell(size(edges,1),1);

if ( isa(bc, 'ultraSEM.BC') )

    if ( length(bc) == 1 )
        % If given a single BC object, replicate it to match the number of
        % boundaries.
        bc = repelem(bc, length(coeffs));
    end

    if ( length(bc) == length(coeffs) )
        % TODO: The use of S.patches{1}.D2N and edges breaks encapsulation.
        D2N = S.patches{1}.D2N;
        nn = sum(edges(:,5));
        A = eye(nn);
        bcval = zeros(nn, 1);
        idx = 0;

        % TODO: Detect pure Neumann boundary conditions.
        for k = 1:size(coeffs, 1)
            if ( isempty(bc(k).dir) || isempty(bc(k).neu) || ...
                 isempty(bc(k).val) || (bc(k).dir == 0 && bc(k).neu == 0) )
                error('ULTRASEM:ULTRASEM:bc2coeffs:unspecified', ...
                    'Unspecified boundary condition at position %i.', k);
            end
            n = edges(k,5);
            side = (idx+1):(idx+n);

            % Construct the operator that maps from Dirichlet data to mixed
            % boundary data on this boundary:
            A(side,:) = bc(k).dir*A(side,:) + bc(k).neu*D2N(side,1:end-1);

            % Get the coefficients of the boundary data:
            cfs = toCoeffs(bc(k).val, edges(k,:));

            % If there was a Neumann contribution on this boundary, we must
            % subtract off the normal derivative of the particular
            % solution.
            bcval(side) = cfs - bc(k).neu*D2N(side,end);
            idx = idx + n;
        end

        % Map the given mixed boundary data to pure Dirichlet data:
        bc = A \ bcval;
    else
        error('ULTRASEM:ULTRASEM:bc2coeffs:length', ...
            ['The number of boundary conditions does not match the' ...
             'number of boundaries.']);
    end

end

if ( ~isnumeric(bc) || isscalar(bc) )
    for k = 1:size(coeffs, 1)
        coeffs{k} = toCoeffs(bc, edges(k,:));
    end
elseif ( isvector(bc) )
    % We were given a vector, so assume these are coefficients.
    idx = 0;
    for k = 1:size(coeffs, 1)
        n = edges(k,5);
        side = (idx+1):(idx+n);
        coeffs{k} = bc(side);
        idx = idx + n;
    end
else
    error('ULTRASEM:ULTRASEM:bc2coeffs:badBC', ...
        'Cannot evaluate boundary data.');
end

% CAT() is 10x faster than CELL2MAT().
coeffs = cat(1, coeffs{:});

end

function coeffs = toCoeffs(bc, edge)
%TOCOEFFS   Convert boundary data to Chebyshev coefficients.
%   TOCOEFFS(BC, EDGE) converts the boundary data BC along the boundary
%   EDGE to a vector of univariate Chebyshev coefficients. BC may be a
%   function handle, a chebfun2, or scalar. EDGE is a row vector that
%   conforms to the structure of ultraSEM.Patch.edges.

if ( ~isnumeric(bc) )
    % Evaluate the BC if given a function handle:
    % Create grid on this edge:
    a = edge(1)+edge(2)*1i;
    b = edge(3)+edge(4)*1i;
    t = chebpts(edge(5));
    x = real(b-a)/2*t+real(b+a)/2;
    y = imag(b-a)/2*t+imag(b+a)/2;
    % Evaluate at grid:
    vals = feval(bc, x, y);
    % Convert from values to coeffs:
    coeffs = util.vals2coeffs(vals);
elseif ( isscalar(bc) )
    % Convert a scalar to coeffs:
    n = edge(5);
    coeffs = [bc ; zeros(n-1, 1)];
end

end
