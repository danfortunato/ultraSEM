function out = quad2rect(varargin)
%QUAD2RECT   Convert from a 4x2 to 1x4 representation of a rectangle.

% Deal with cell input:
if ( nargin == 1 )
    v = varargin{1};
else
    v = varargin;
end
if ( iscell(v) )
    nv = numel(v);
    out = zeros(nv,4);
    for k = 1:nv
        out(k,:) = util.quad2rect(v{k});
    end
    return
end       

% Find max in x and y
x1 = min(v(:,1));
x2 = max(v(:,1));
y1 = min(v(:,2));
y2 = max(v(:,2));
% 4x1 representation:
out = [x1, x2, y1, y2];

end
