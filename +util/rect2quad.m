function out = rect2quad(v)
%RECT2QUAD   Convert from a 1x4 to 4x2 representation of a rectangle.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

if ( isempty(v) )
    out = [];
    return
end

% Deal with cell input:
if ( size(v, 2) ~= 4 )
    v = v.';
end

% v may be a nx4 matrix of rectangles:
out = {};
for k = 1:size(v, 1)
    vk = v(k,:);
    out{k} = vk([1 3 ; 2 3 ; 2 4 ; 1 4]);
end

if ( numel(out) == 1 )
    out = out{:};
end

end