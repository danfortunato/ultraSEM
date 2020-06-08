function n = normals(Q)
%NORMALS   Outward pointing normal vectors to the edges of a mapping.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

v = [ Q.x([-1 1 1 -1],[-1 -1 1 1]) ;
      Q.y([-1 1 1 -1],[-1 -1 1 1]) ].';
v = [ v ; v(1,:) ];
dx = diff(v(:,1));
dy = diff(v(:,2));
n = [dy.' ; -dx.' ];
n = n(:,[4 2 1 3]);  % Reorder to left, right, down, up
n = n * diag(1./sqrt(sum( n.^2 )));  % Normalize

end
