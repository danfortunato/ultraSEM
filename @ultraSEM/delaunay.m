function T = delaunay(v, P)
%DELAUNAY   Delaunay triangulation.
%   T = ULTRASEM.DELAUNAY(X), where is an N x 2 matrix, returns an
%   ULTRASEM.DOMAIN representation of the Delaunay triangulation equivalent
%   to DELAUNAY(X), but with each triangle subdivided into three further
%   quadrilaterals with their intersection at the barycenter.
%
%   T = ULTRASEM.DELAUNAY(X_BDY, X_INT), where X_INT is an N x 2 matrix and
%   X_BDY is an M x 2 matrix, creates a constrained Delaunay triangulation,
%   equivalent to delaunayTriangulation(X_BDY, X_INT).
%
% See also TRIANGLE, QUAD.

% Compute triangulation:
if ( nargin < 2 )
    dt = delaunayTriangulation(v);
    list = dt.ConnectivityList;
else
    dt = delaunayTriangulation(v, P);
    IO = isInterior(dt);
    list = dt(IO,:);
end

pts = dt.Points(:,:);
% Build domain:
nt = size(list,1);
T = [];
for k = 1:nt
    v = pts(list(k,:),:);
    T = T & ultraSEM.triangle(v);
end

end
