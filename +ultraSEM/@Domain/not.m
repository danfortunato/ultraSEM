function E = not(T, pad)
%~   Take the complement of an ULTRASEM.DOMAIN.
%   ~T forms a rectangular domain by constructing a rectangular domain R
%   containing T and removing T from it. The size of R will be one patch
%   size greater than the extremities of T.
%
%   NOT(T, P) is similar to the above, but will pad by an amount P on each
%   side of T. If P is a 4x1 vector then it is interpreted as
%   P = [P_left, P_right, P_bottom, P_top].

% Determine padding size:
if ( nargin == 1 )
    pad = [1 1 1 1];
elseif ( isscalar(pad) )
    pad = repmat(pad, 1, 4);
end

% Determine sizes:
dom = T.domain;
dx = diff(dom(1,1:2)); dy = diff(dom(1,3:4));
x1 = min(dom(:,1));    x2 = max(dom(:,2));
y1 = min(dom(:,3));    y2 = max(dom(:,4));

% Construct containing rectangle:
newDom = [x1-pad(1)*dx, x2+pad(2)*dx y1-pad(3)*dy y2+pad(4)*dy];
m = diff(newDom(1:2))/dx; n = diff(newDom(3:4))/dy;
R = ultraSEM.rectangle(newDom, m, n);

% Extract T from R:
E = R - T;

end
