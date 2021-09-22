function S = reorderEdges(S)

if ( numel(S.patches) > 1)
    build(S);
end

edges = S.patches{1}.edges;
n = size(edges,1);

initialOrdering = 1:n; 
p = [1 ; zeros(n-1,1)];
unassignedEdges = [false ; true(n-1,1)];

for k = 2:n
    [~,loc_k] = ismember(edges(p(k-1),[3 4]), ...
        edges(unassignedEdges,[1 2]), 'rows');
    if ( ~loc_k )
        [~,loc_k] = ismember(edges(p(k-1),[3 4]), ...
            edges(unassignedEdges,[3 4]), 'rows');
    end
    if ( ~loc_k )
        [~,loc_k] = ismember(edges(p(k-1),[1 2]), ...
            edges(unassignedEdges,[3 4]), 'rows');
    end    
    if ( ~loc_k )
        error('uh oh');
    end 
    edgeList = initialOrdering(unassignedEdges);
    p(k) = edgeList(loc_k);
    unassignedEdges(edgeList(loc_k)) = false;
end

<<<<<<< Updated upstream
<<<<<<< Updated upstream
p
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
S.patches{1}.edges = edges(p,:);

end


