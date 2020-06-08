function v = assertIsRect( v )
%ASSERTISRECT   Check we have valid vertices for an ULTRASEM.RECT.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

if ( ~isnumeric(v) )
    error('Input should be numeric.');
end
% Ensure v is of the form [x, y], not [x ; y]:
if ( size(v,2) ~= 2 ), v = v.'; end
% Check dimension:
if ( size(v,2) ~= 2 || size(v, 1) ~= 4 )
    error('Incorrect vertices dimension.')
end
% Check it's a valid rectangle:
ndx = ~diff([v(:,1) ; v(1,1)]);
ndy = ~diff([v(:,2) ; v(1,2)]);
s1 = [0;1;0;1];
s2 = [1;0;1;0];
if ( ~( all(ndx==s1 & ndy==s2) || all(ndx==s2 & ndy==s1) ) )
    error('Not a valid rectangle.')
end

end
