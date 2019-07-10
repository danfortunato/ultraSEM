function T = trimesh(p, t)
%ULTRASEM.TRIMESH   Return a triangulated mesh.
%   T = ULTRASEM.TRIMESH(TRI) returns an ultraSEMDomain mesh T made from
%   the triangulation TRI.
%
%   T = ULTRASEM.TRIMESH(P, T) returns an ultraSEMDomain mesh T made up of
%   triangles specified by vertices P and connectivity list T.

    if ( nargin == 1 && isa(p, 'triangulation') )
        t = p.ConnectivityList;
        p = p.Points;
    end

    if ( size(p,2)~=2 || size(t,2)~=3 )
        error('ULTRASEM:ULTRASEMDOMAIN:trimesh:invalid', ...
            'Invalid triangular mesh specification.');
    end

    T = ultraSEM.triangle(p(t(1,:),:));
    for k = 2:size(t,1)
        T = T & ultraSEM.triangle(p(t(k,:),:));
    end

end
