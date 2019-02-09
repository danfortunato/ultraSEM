function varargout = gui()
% Solve PDO in a domain drawn on a figure.

if ( nargin == 0 )
    PDO = {{1, 0, 1}, {0, 0}, 0};
end
if ( nargin < 2 )
    n = 21;
end
rhs = -1;
bc = 0;

% TODO: Once ultraSEM is more general this can return the domain without
% knowing the PDO in advance.

close all;
LW = 'linewidth';  lw = 2;
MS = 'markersize'; ms = 20;

%% Acquire boundary vertices:

% Launch a new figure and click the mouse as many time as you want within
% figure
lim = 5*[-1 1 -1 1];
figure(), axis(lim), grid on, hold on
title('Select boundary vertices')
v = [];
cv = ginput(1);                             % Current vertex.
v = cv;                                     % Add to list.
plot(cv(1), cv(2), 'b.', MS, ms), drawnow   % Plot vertex
button = 1;
while ( button == 1  )                      % Until button 2 is pressed...
    cv = v(end,:);                          % Current vertex.
    [x, y, button] = ginput(1);             % Grab next vertex.
    if ( button == 1 )
        nv = [x, y];                        % Next vertex.
    else
        nv = v(1, :);                       % Close the domain.
    end
    plot(nv(1), nv(2), 'b.', MS, ms)        % Plot next vertex.
    ce = [cv ; nv];                         % Current edge.
    plot(ce(:, 1), ce(:, 2), 'b-', LW, lw)  % Draw current edge.
    v = [v ; nv];                           % Add to list.
    if all( nv == v(1,:) )                  % Break when domain is closed.
        break
    end
end

% Check that the domain is closed.
if ( ~all( v(end,:) == v(1,:) ) )
    error('ultraSEMGUI:OpenDomain', 'Domain must be closed.')
end
% Remove the duplicated end vertex:
v(end, :) = [];

%% Determine polygonal boundary contraints:

lv = size(v, 1);
P = [1:lv ; 2:lv 1].';

%% Add interior vertices:

title('Select interior vertices')
button = 1;
while ( button == 1  )                  % Until button 2 is pressed...
    [x, y, button] = ginput(1);         % Grab next vertex.
    if ( button == 1 )
        nv = [x, y];                    % Next vertex.
    else
        break
    end
    plot(nv(1), nv(2), 'r.', MS, ms)    % Plot next vertex.
    v = [v ; nv];                       % Add to list.
end

%% Compute and plot triangularization:

title('Triangularization')
dt = delaunayTriangulation(v, P);
IO = isInterior(dt);
list = dt(IO,:);
pts = dt.Points(:,:);
triplot(list, pts(:,1), pts(:,2))
drawnow, hold off

%% Build domain:

nt = size(list,1);
T = [];
for k = 1:nt
    v = pts(list(k,:),:);
    T = T & ultraSEM.triangle(v);
end

%% Initialize ultraSEM object:
S = ultraSEM(T, PDO, rhs, n);

%% Solve:
sol = solve(S, bc);
if ( nargout == 1 )
    varargout{1} = sol;
end

%% Plot:
title('Solution')
figure(1)
h = surf(sol);
axis(lim), axis(h(1).Parent, lim)

end
