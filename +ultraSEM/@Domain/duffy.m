function T = duffy(varargin)

if ( nargin == 0 )
    vertices = [0 0; 1 0; 0 1];
else
    vertices = varargin{1};

    if ( ultraSEM.Domain.isClockwise(vertices) )
        % Switch the second and third indices.
        vertices([2,3],:) = vertices([3,2],:);
    end
end

K = ultraSEM.Tri( vertices );
T = ultraSEM.Domain(K);

end
