function T = polygon(v)
%ULTRASEM.POLYGON  Return a convex polygonal domain formed of kites.
%   T = ULTRASEM.POLYGON(V) returns a convex polygonal ultraSEMDomain T
%   with vertices V formed of N kites, where N is the polygon degree. The
%   vertices V should be given in anticlockwise order. If they are not,
%   they will be modified to be as such.
%
%   It is important to note that each side of the polygon will be formed
%   from two adjacent grids.

    % Number of sides of the polygon
    n = size(v,1);

    % Check that the given points form a convex polygon:
    if ( ~all(unique(convhull(v)) == (1:n)') )
        error('ULTRASEM:ULTRASEMDOMAIN:polygon:nonconvex', ...
            'Polygon is not convex.');
    end

    if ( ultraSEMDomain.isClockwise(v) )
        % Switch the vertices to make it anticlockwise.
        v = flipud(v);
    end

    % Locate the center of the polygon:
    c = mean( v );

    % Construct kites:
    K(n,1) = kite();
    for k = 1:n
        next = mod(k,n)+1;
        prev = mod(k-2,n)+1;
        K(k) = kite( [ v(k,:) ; mean(v([k,next],:)) ; c ; mean(v([prev,k],:)) ] );
    end

    % Construct merges:
    lvls = ceil(log2(n));
    mergeIdx = cell(1, lvls);
    for k = 1:lvls
        m = ceil(2^(1-k)*n);
        mergeIdx{k} = reshape( padarray((1:m)', mod(m,2), NaN, 'post'), 2, ceil(m/2))';
    end

    % Construct ultraSEMDomain:
    T = ultraSEMDomain(K, mergeIdx);

end
