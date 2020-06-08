function plottree(T)

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

%TODO: Experimental
dom = T.domain;
idx = T.mergeIdx;

for k = 1:numel(dom)
    newDom{k} = dom(k,:);
end
for k = 1:min(numel(idx), 9)
    figure(k)
    dom = newDom;
    idxk = idx{k};
    for j = 1:size(idxk,1)
        if ( isnan(idxk(j,2)) )
            newDom{j} = dom{idxk(j,1)};
        else
            newDom{j} = [dom{idxk(j,1)}, dom{idxk(j,2)}];
        end
        plot(newDom{j}); hold on
    end
    hold off
    alignfigs
end

end
