function v = assertIsTri(v)
%ASSERTISTRI   Check we have valid vertices for an ULTRASEM.TRI.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

if ( ~isnumeric(v) )
    error('Input should be numeric.');
end
% Ensure v is of the form [x, y], not [x ; y]:
if ( size(v,2) ~= 2 ), v = v.'; end
% Check dimension:
if ( size(v,2) ~= 2 || size(v, 1) ~= 3 )
    error('Incorrect vertices dimension.')
end
% Ensure vertices are oriented in an anticlockwise direction:
if ( ultraSEM.Domain.isClockwise(v) )
    % Switch the second and third indices.
    v([2,3],:) = v([3,2],:);
end

end
