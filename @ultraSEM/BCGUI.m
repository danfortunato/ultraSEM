function bc = BCGUI(S, h)

bc = [];
if ( nargin == 2 )
    % h = h;
elseif ( isa(S, 'ultraSEM') )
    h = plot(S); shg
else
    % This should not really be supported..
    h = S;
end
set(gca, 'ButtonDownFcn', @myCallback);
set(gca, 'DeleteFcn', @deleteCallback);
drawnow
uiwait(gcf)
    
    function deleteCallback(varargin)
        bc = ultraSEM.BC.getBCsFromPlot(h);
    end

    function myCallback(varargin)
        
        % Select bounding box:
        a = gca;
        point1 = a.CurrentPoint;           % button down detected
        rbbox;                             % initialise rubberband box
        point2 = a.CurrentPoint;           % button up detected
        point1 = point1(1,1:2);            % extract x and y
        point2 = point2(1,1:2);
        p1 = min(point1,point2);           % calculate locations
        offset = abs(point1-point2);       % and dimensions
        vx = [p1(1) p1(1)+offset(1) p1(1)+offset(1) p1(1) p1(1)];
        vy = [p1(2) p1(2) p1(2)+offset(2) p1(2)+offset(2) p1(2)];

        % Locate chosen edges:
        xy = vertcat(h.Position);
        x = xy(:,1); y = xy(:,2);
%         idx = inpolygon(x, y, vx, vy);
        idx = insquare(x, y, vx, vy);

        if ( any(idx) )
            idx = find(idx);
            % Aquire values:
            data = inputdlg({'Dirichlet','Neumann','Value'},...
                ['Boundary condition ' h(idx(1)).String], ...
                [1 50], h(idx(1)).UserData);
            if ( isempty(data) )
                return
            end
            for k = 1:numel(idx)
                % Set data
                h(idx(k)).UserData(1:3) = data;
                % Adjust linestyle:
                if (data{1} == '0')
                    set(h(idx(k)).UserData{4}, 'LineStyle', ':');
                else
                    set(h(idx(k)).UserData{4}, 'LineStyle', '-');
                end
            end
        end

    end

    function idx = insquare(x,y,vx,vy)

        idxx = (x > min(vx)) & (x < max(vx));
        idxy = (y > min(vy)) & (y < max(vy));
        idx = idxx & idxy;

    end

end