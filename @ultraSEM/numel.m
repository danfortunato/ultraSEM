function N = numel(S)
%NUMEL   Number of degrees of freedom in an ULTRASEM.
%   N = NUMEL(S) returns the total number of degrees of freedom in the
%   ULTRASEM object S. If S is an array of objects, then N is the number of
%   objects in the array.

numObj = builtin('numel', S);
if ( numObj > 1 )
    % This is an array of objects. Don't overload.
    N = numObj;
    return
end

N = 0;
for k = 1:numel(S.patches)
    N = N + numel(S.patches{k});
end

end
