function out = rect2quad(v)

if ( size(v, 2) ~= 4 )
    v = v.';
end


out = {};
for k = 1:size(v, 1)
    vk = v(k,:);
    out{k} = vk([1 3 ; 2 3 ; 2 4 ; 1 4]);
end

if ( numel(out) == 1 )
    out = out{:};
end

end