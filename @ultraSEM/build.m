function build(S)
%BUILD   Build the merge tree of patches in an ULTRASEM object.
%   BUILD(S) will construct the global solution operator from the local
%   operators in S.patches in the order defined in the merge tree
%   S.domain.mergeIdx. After building, S will have a single patch.
%
%   If S has not yet been initialized, then an error is thrown.
%
% See also INITIALIZE, SOLVE.

if ( isa(S.patches, 'ultraSEM') )
    % Recurse down and build lower level patches
    for k = 1:numel(S.patches)
        build(S.patches(k));
    end
    % Concatenate (since S already contains domain info.)
    S.patches = vertcat(S.patches.patches); % TODO: breaks encapsulation
end

if ( ~isInitialized(S) )
    error('ULTRASEM:ULTRASEM:build:notInitialized', ...
        '%f has not yet been initialized.', inputname(1))
end

% Build the patches:
S.patches = build(S.domain, S.patches);

end
