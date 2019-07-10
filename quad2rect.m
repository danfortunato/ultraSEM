function out = quad2rect(v)

x1 = min(v(:,1));
x2 = max(v(:,1));
y1 = min(v(:,2));
y2 = max(v(:,2));

out = [x1, x2, y1, y2];

end