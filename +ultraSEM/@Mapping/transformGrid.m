function [X, Y, XY] = transformGrid(T, x, y)
%TRANSFORMGRID   Map a grid.

X = T.x(x, y);
Y = T.y(x, y);

if ( nargout == 3 )
    XY = [X Y];
end

end
