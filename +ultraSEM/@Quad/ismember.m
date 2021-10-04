function idx = ismember(z, Q)
%ISMEMBER  Check if points are within a Quad.
%   ISMEMBER(Z, Q) returns true if the point Z=[X,Y] is (strictly) within
%   the quad Q and false otherwise. Z may may be an Nx2 matrix, in which
%   case ISMEMBER returns an Nx1 vector True/False values.

% TODO: It's maybe faster to use the inverse mapping to see if m^{-1}(z)
% lies within the reference square.

v = Q.v;
vx = v(:,1);
vy = v(:,2);
x = z(:,1);
y = z(:,2);

[idx, on] = inpolygon(x, y, vx, vy);
idx(on) = false;

end