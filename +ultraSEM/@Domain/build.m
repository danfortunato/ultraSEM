function p = build(T, p)

% TODO: Document.

if ( length(T) == 1 || numel(p) == 1 )
    % Nothing to build for a single patch.
    return
elseif ( isempty(T.mergeIdx) )
    warning('ULTRASEM:DOMAIN:build:noMergeIdx', ...
        '%s does not contain merge information.', inputname(1))
    T.mergeIdx = T.defaultIdx(T.domain);
end

% Loop over each level:
for j = 1:numel(T.mergeIdx) % <-- number of levels

    idxj = T.mergeIdx{j};         % Indicies for jth level.
    q = cell(size(idxj, 1),1);    % Initialise patches at jth level.

    % Perform all the merges at this level:
    for k = 1:size(idxj, 1)
        idxjk = idxj(k,:);        % Index for kth merge at jth level.
        idxjk(isnan(idxjk)) = []; % Remove NaNs.
        pk = p(idxjk);            % Extract patches.
        pk(cellfun(@isempty, pk)) = []; % Remove empty patches.
        if ( isempty(pk) )
            q{k} = [];            % Empty merge.
        else
            q{k} = merge(pk{:});  % Perform merge.
        end
    end

    p = q; % Store back in p for recursion.

end

end
