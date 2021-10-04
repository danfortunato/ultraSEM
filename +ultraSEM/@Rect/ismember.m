function idx = ismember(z, Q)
%ISMEMBER  Check if points are within a Rect.
%   ISMEMBER(Z, Q) returns true if the point Z=[X,Y] is (strictly) within
%   the Rect Q and false otherwise. Z may may be an Nx2 matrix, in which
%   case ISMEMBER returns an Nx1 vector True/False values.

v = Q.v;
vx = v(:,1);
vy = v(:,2);
x = z(:,1);
y = z(:,2);

idxx = (x > min(vx)) & (x < max(vx));
idxy = (y > min(vy)) & (y < max(vy));
idx = idxx & idxy;

end