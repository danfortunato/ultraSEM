function out = quad2rect(varargin)

if ( nargin == 1 )
    v = varargin{1};
else
    v = varargin;
end

if ( iscell(v) )
    nv = numel(v);
    out = zeros(nv,4);
    for k = 1:nv
        out(k,:) = quad2rect(v{k});
    end
    return
end       

x1 = min(v(:,1));
x2 = max(v(:,1));
y1 = min(v(:,2));
y2 = max(v(:,2));

out = [x1, x2, y1, y2];

end