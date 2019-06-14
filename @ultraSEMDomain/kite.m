function T = kite(varargin)

vertices = varargin{1};

if ( ultraSEMDomain.isClockwise(vertices) )
    % Switch the second and fourth indices.
    vertices([2,4],:) = vertices([4,2],:);
end

K = kite( vertices );
%T = ultraSEMDomain(K, {[1 ; NaN]});
T = ultraSEMDomain(K);

end
