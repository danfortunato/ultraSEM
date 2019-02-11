function T = pentaflake(n, varargin)
%PENTAFLAKE   Create a Penrose snowflake domain, iterated N times.

    optargs = {0 0 6*pi/5 1};
    optargs(1:nargin-1) = varargin;
    [x, y, th, h] = optargs{:};

    phi = (1+sqrt(5))/2;
    tt = linspace(th, th+2*pi, 6).'; tt(end) = [];
    [offx, offy] = pol2cart(tt+pi/10, h*phi);

    if ( n == 0 )
        T = ultraSEM.polygon([h*sin(tt)+x, h*cos(tt)+y]);
    else
        if ( n > 1 ), h = h/(phi+1); end
        T = [];
        for k = 1:5
            T = T & pentaflake(n-1, x+offx(k), y+offy(k), th, h);
        end
        T = T & pentaflake(n-1, x, y, th+pi, h);
    end

end
