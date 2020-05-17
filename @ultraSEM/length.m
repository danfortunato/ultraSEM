function n = length(S)
%LENGTH   Number of patches in an ULTRASEM.
%   LENGTH(S) returns the total number of patches in the ULTRASEM object S.
%   If S is an array of objects, then LENGTH(S) is the number of objects in
%   the array.

numObj = builtin('length', S);
if ( numObj > 1 )
    % This is an array of objects. Don't overload.
    n = numObj;
    return
end

n = 0;
for k = 1:numel(S.patches)
    n = n + length(S.patches{k});
end

end
