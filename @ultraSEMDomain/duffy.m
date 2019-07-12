function T = duffy(varargin)

if ( nargin == 0 )
    vertices = [];
else
    vertices = varargin{1};

    if ( ultraSEMDomain.isClockwise(vertices) )
        % Switch the second and third indices.
        vertices([2,3],:) = vertices([3,2],:);
    end
end

K = ultraSEMTri( vertices );
T = ultraSEMDomain(K);

end
