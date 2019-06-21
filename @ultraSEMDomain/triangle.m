function T = triangle(v)
%ULTRASEM.TRIANGLE  Return a triangular domain formed of kites/quads.
%   T = ULTRASEM.TRIANGLE(V) returns a triangular ultraSEMDomain T with
%   vertices V formed of three kites. The vertices V should be given in
%   anticlockwise order. If they are not, they will be modified to be as
%   such.
%
%   It is important to note that each side of the triangle will be formed
%   from two adjacent grids.

    if ( nargin == 0 )
        v = [0 0 ; 1 0 ; .5 sqrt(3)/2];
    end

    if ( ultraSEMDomain.isClockwise(v) )
        % Switch the second and third vertex to make it anticlockwise.
        v([2,3],:) = v([3,2],:);
    end

    % Locate the centre of the triangle:
    c = mean( v );

    % Determine vertices of kites:
    m12 = mean( v([1,2],:) );
    m13 = mean( v([1,3],:) );
    m23 = mean( v([2,3],:) );
    v1 = [ v(1,:) ; m12 ; c ; m13];
    v2 = [ m13 ; c ; m23 ; v(3,:) ];
    v3 = [ m12 ; v(2,:) ; m23 ; c ];

    % Construct kites:
    K(3,1) = ultraSEMQuad();
    K(1) = ultraSEMQuad( v1 );
    K(2) = ultraSEMQuad( v2 );
    K(3) = ultraSEMQuad( v3 );

    % Construct ultraSEMDomain:
    T = ultraSEMDomain(K, {[1 2 ; 3 NaN], [1 2]});

end
