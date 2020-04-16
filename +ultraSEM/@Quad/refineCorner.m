function [Q, idx] = refineCorner(Q, k)
%REFINECORNER   Refine an ULTRASEM.QUAD into a corner.
%
% See also REFINE, REFINEPOINT.

% TODO: Transpose everythng.
v = vertices(Q)';
c = centroid(Q);
Q(3,1) = ultraSEM.Quad();
if ( k == 1 )
    vnew = (v(:,1) + v(:,[2,4]))/2;
    v1 = [v(:,1) vnew(:,1) c vnew(:,2)];
    v2 = [v(:,[2 3]) c vnew(:,1)];
    v3 = [v(:,4) vnew(:,2) c v(:,3)];
elseif ( k == 2 )
    vnew = (v(:,2) + v(:,[3,1]))/2;
    v1 = [v(:,2) vnew(:,1) c vnew(:,2)];
    v2 = [v(:,[3 4]) c vnew(:,1)];
    v3 = [v(:,1) vnew(:,2) c v(:,4)];                
elseif ( k == 3 )
    vnew = (v(:,3) + v(:,[4,2]))/2;
    v1 = [v(:,3) vnew(:,1) c vnew(:,2)];
    v2 = [v(:,[4 1]) c vnew(:,1)];
    v3 = [v(:,2) vnew(:,2) c v(:,1)]; 
elseif ( k == 4 )
    vnew = (v(:,4) + v(:,[1,3]))/2;
    v1 = [v(:,4) vnew(:,1) c vnew(:,2)];
    v2 = [v(:,[1 2]) c vnew(:,1)];
    v3 = [v(:,3) vnew(:,2) c v(:,2)];                 
end 
Q(3) = ultraSEM.Quad(v1');
Q(2) = ultraSEM.Quad(v2');            
Q(1) = ultraSEM.Quad(v3');
idx = {[1 2 ; 3 NaN], [1 2]};

end

% function Q = refineCorner(Q, k)
%     if ( k ~= 1)
%         error ('only k = 1 implemented')
%     end
%     v = vertices(Q);
%     c = centroid(Q).';
%     vnew = (v(1:end,:) + v([2:end, 1],:))/2;
%     Q(3,1) = ultraSEM.Quad();
%     Q(1) = ultraSEM.Quad([v(1,:) ; vnew(1,:) ; c ; vnew(4,:)]);
%     Q(2) = ultraSEM.Quad([v(2,:) ; v(3,:) ; c ; vnew(1,:)]);            
%     Q(3) = ultraSEM.Quad([v(4,:) ; vnew(4,:) ; c ; v(3,:)]);
% end
