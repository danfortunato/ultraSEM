function [Q, idx] = refinePoint(Q, z)
%REFINEPOINT   Refine an ULTRASEM.QUAD around a point.
%
%   See also REFINE, REFINECORNER.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

idx = {[1, NaN]};

% Check to sese if the point is a corner:
loc_corner = find(ismember(Q.v, z, 'rows'));
% Check to see if the point is interior:
loc_interior = find(ismember(z, Q));

% TODO: Need to fix this...
if ( (numel(loc_corner) > 1) || (numel(loc_interior) > 1) ...
        || (any(loc_interior) && any(loc_corner)) )
    error('Can only refine around a single point per patch.');
end

if ( any(loc_corner) )
    % Refine a corner point:
    [Q, idx] = refineCorner(Q, loc_corner);
    
elseif ( any(loc_interior) )
    % REfine an interior point:
        
    z = z(loc_interior,:)';
    v = vertices(Q)';
    vnew = (z + v)/2;
    v1 = [v(:,[1 2]) vnew(:,[2 1])];
    v2 = [v(:,[2 3]) vnew(:,[3,2])];
    v3 = [v(:,[3 4]) vnew(:,[4 3])];
    v4 = [v(:,[4 1]) vnew(:,[1 4])];
    
    Q(5,1) = ultraSEM.Quad();
    Q(1) = ultraSEM.Quad(v1');
    Q(2) = ultraSEM.Quad(v2');
    Q(3) = ultraSEM.Quad(v3');            
    Q(4) = ultraSEM.Quad(v4');            
    Q(5) = ultraSEM.Quad(vnew');   

    idx = {[1 2 ; 3 4 ; 5 NaN], [1 2 ; 3 NaN] , [1 2]};    
end

   
end