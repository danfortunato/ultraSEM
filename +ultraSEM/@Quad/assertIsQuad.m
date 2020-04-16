function v = assertIsQuad(v)
%ASSERTISQUAD   Check we have valid vertices for an ULTRASEM.QUAD.

if ( ~isnumeric(v) )
    error('Input should be numeric.');
end
% Ensure v is of the form [x, y], not [x ; y]:
if ( size(v,2) ~= 2 ), v = v.'; end
% Check dimension:
if ( size(v,2) ~= 2 || size(v, 1) ~= 4 )
    error('Incorrect vertices dimension.')
end
% Ensure vertices are oriented in an anticlockwise direction:
if ( ultraSEM.Domain.isClockwise(v) )
    % Switch the second and fourth indices.
    v([2,4],:) = v([4,2],:);
end

end
