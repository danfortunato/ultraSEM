function T = rectangle(dom, m, n)
%RECTANGLE   Return a rectangular domain.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

if ( nargin == 2 && strcmp(m, 'makeObj') )
    n = size(dom, 1);
    T(n,1) = rectangle(dom(n,:));
    for k = 1:n-1
        T(k) = rectangle(dom(k,:));
    end
    return
end

% Default rectangle
if ( nargin == 0 )
    dom = [-1 1 -1 1];
end

if ( ~all(size(dom) == [1, 4]) )
    error('ULTRASEM:DOMAIN:rectangle', ...
        'Input must be a 1x4 vector.')
end
if ( dom(1) == dom(2) || dom(3) == dom(4) )
    error('ULTRASEM:DOMAIN:rectangle', ...
        'Rectangle is degenerate.');
end
if ( nargin < 2 )
    m = 1;
end
if ( nargin < 3 )
    n = m;
end

% Initialize tree:
idx = {};

if ( m == 1 && n == 1 )
    % Only a single patch. No subdivision or merge info required.
    dom = ultraSEM.Rect(dom);
    T = ultraSEM.Domain(dom, {[1 NaN]});
    return
end

x = linspace(dom(1), dom(2), n+1);
y = linspace(dom(3), dom(4), m+1);
l = 0;
for j = 1:m
    for k = 1:n
        l = l+1;
        dom(l,:) = [x(k), x(k+1), y(j), y(j+1)];
    end
end

if ( (log2(n) == round(log2(n))) && (log2(m) == round(log2(m))) )
    % Powers of two are easier to contend with

    l = 1;
    while ( n*m > 1 )
        tmp = reshape(1:n*m, n, m);
        if ( n >= m )
            % Merge up
            n = n/2;
        else
            % Merge right
            tmp = tmp.';
            m = m/2;
        end
        tmp = reshape(tmp, 2, n*m)';
        idx{l} = tmp; l = l + 1;
    end

else

    l = 0;
    while ( n*m > 1 )
        l = l+1;
        idx{l} = [];
        if ( m <= n )   % Merge right
            for k = 1:m
                idxk = (k-1)*n + (1:n);
                if ( mod(n, 2) )
                    idxk = [idxk, NaN]; %#ok<AGROW>
                end
                idxk = reshape(idxk, 2, length(idxk)/2)';
                idx{l} = [idx{l} ; idxk];
            end
            n = ceil(n/2);
        else            % Merge up
            for k = 1:ceil(m/2)
                idxk = 2*(k-1)*n + (1:2*n);
                idxk(idxk > n*m) = NaN;
                idxk = reshape(idxk, length(idxk)/2, 2);
                idxk(all(isnan(idxk), 2),:) = [];
                idx{l} = [idx{l} ; idxk];
            end
            m = ceil(m/2);
        end
    end

end

dom = ultraSEM.Rect(dom);
T = ultraSEM.Domain(dom, idx);

end
