function T = flatten(T)
%FLATTEN   Try to flatten merges (for efficiency).

if ( numel(T) > 1 )
    for k = 1:numel(T)
        T(k) = flatten(T(k));
    end
    flattenthis_TODO
    return
end
for k = 1:numel(T.domain)
    if ( isa(T.domain(k), 'ultraSEM.Domain') && ...
            isa(T.domain(k).domain, 'ultraSEM.Domain') )
        T.domain(k) = flatten(T.domain(k));
    end
end

if ( isRect(T.domain(1).domain(1)) == isRect(T.domain(2).domain(1)) )
    T = oldMerge(T.domain(1), T.domain(2));
end        

end
