function T = quad(varargin)

vertices = varargin{1};

if ( ultraSEM.Domain.isClockwise(vertices) )
    % Switch the second and fourth indices.
    vertices([2,4],:) = vertices([4,2],:);
end

K = ultraSEM.Quad( vertices );
%T = ultraSEM.Domain(K, {[1 ; NaN]});
T = ultraSEM.Domain(K);

end
