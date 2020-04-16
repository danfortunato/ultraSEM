function T = triangle(v)
%ULTRASEM.TRIANGLE   Return a triangular domain formed of quadrilaterals.
%   T = ULTRASEM.TRIANGLE(V) returns a triangular ULTRASEM.DOMAIN T with
%   vertices V formed of three kites. The vertices V should be given in
%   anticlockwise order. If they are not, they will be modified to be as
%   such.
%
%   It is important to note that each side of the triangle will be formed
%   from two adjacent grids.

    if ( nargin == 0 )
        v = [0 0 ; 1 0 ; .5 sqrt(3)/2];
    end

    if ( ultraSEM.Domain.isClockwise(v) )
        % Switch the second and third vertex to make it anticlockwise.
        v([2,3],:) = v([3,2],:);
    end

    % Locate the centre of the triangle:
    c = mean( v );

    % Determine vertices of quads:
    m12 = mean( v([1,2],:) );
    m13 = mean( v([1,3],:) );
    m23 = mean( v([2,3],:) );
    v1 = [ v(1,:) ; m12 ; c ; m13];
    v2 = [ m12 ; v(2,:) ; m23 ; c ];
    v3 = [ m13 ; c ; m23 ; v(3,:) ];

    % Construct quads:
    K(3,1) = ultraSEM.Quad();
    K(1) = ultraSEM.Quad( v1 );
    K(2) = ultraSEM.Quad( v2 );
    K(3) = ultraSEM.Quad( v3 );

    % Construct domain:
    T = ultraSEM.Domain(K, {[1 2 ; 3 NaN], [1 2]});

end
