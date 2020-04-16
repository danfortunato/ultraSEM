function [X, Y, XY] = transformGrid(T, x, y)
%TRANSFORMGRID   Map points in [-1 1]^2 to the domain of an ULTRASEM.RECT.

rect = rectVertices(T);
domx = rect(1:2);    domy = rect(3:4); 
sclx = 2/diff(domx); scly = 2/diff(domy);
X = (x+1)/sclx + domx(1); 
Y = (y+1)/scly + domy(1);
if ( nargout == 3 ), XY = [X Y]; end

end
