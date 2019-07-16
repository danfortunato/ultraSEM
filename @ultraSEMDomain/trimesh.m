function T = trimesh(p, t)
%ULTRASEM.TRIMESH   Return a triangulated mesh.
%   T = ULTRASEM.TRIMESH(TRI) returns an ultraSEMDomain mesh T made from
%   the triangulation TRI.
%
%   T = ULTRASEM.TRIMESH(P, T) returns an ultraSEMDomain mesh T made up of
%   triangles specified by vertices P and connectivity list T.

    if ( nargin == 1 && isa(p, 'triangulation') )
        tri = p;
        t = p.ConnectivityList;
        p = p.Points;
    else
        tri = triangulation();
        tri.ConnectivityList = t;
        tri.Points = p;
    end

    if ( size(p,2)~=2 || size(t,2)~=3 )
        error('ULTRASEM:ULTRASEMDOMAIN:trimesh:invalid', ...
            'Invalid triangular mesh specification.');
    end

    T = ultraSEM.triangle(p(t(1,:),:));
    for k = 2:size(t,1)
        T = T & ultraSEM.triangle(p(t(k,:),:));
    end
    
    T.mergeIdx(3:end) = [];
    A = generateAdjacencyMatrix(tri);
    idx = generateMergeIdx(A);
    T.mergeIdx = [T.mergeIdx , idx];
    
end

function A = generateAdjacencyMatrix(tri)
numtri = size(tri,1);
n = neighbors(tri);
A = sparse(numtri,numtri);
for j = 1:numtri
    nk = n(j,:);
    nk(isnan(nk)) = [];
    A(j,nk) = 1;
end

end

function idx = generateMergeIdx(A)

idx = {};
while ( numel(A) > 1 )
    jk = [];
    Anew = A;
    for j = 1:size(A, 1)
        if ( ~isempty(jk) && ismember(j, jk(:,2)) )
            continue
        end
        k = find(A(j,:), 1, 'first');
        if ( isempty(k) )
            jk = [jk ; j NaN];
        else
            jk = [jk ; j k];
            A(:, [j k]) = 0;
            A([k j], :) = 0;
            Anew(j,:) = Anew(j,:) | Anew(k,:);
            Anew(:,j) = Anew(:,j) | Anew(:,k);
            Anew(j,j) = 0;
        end
    end
    idx = [idx , jk];
    A = Anew(jk(:,1),jk(:,1));
end
end

