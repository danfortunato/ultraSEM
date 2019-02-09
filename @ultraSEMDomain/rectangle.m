function T = rectangle(dom, m, n)
%RECTANGLE   Return a rectangular domain.

if ( nargin == 2 && strcmp(m, 'makeObj') )
    n = size(dom, 1);
    T(n,1) = rectangle(dom(n,:));
    for k = 1:n-1
        T(k) = rectangle(dom(k,:));
    end
    return
end

if ( ~all(size(dom) == [1, 4]) )
    error('ULTRASEM:ULTRASEMDOMAIN:rectangle', ...
        'Input must be a 1x4 vector.')
end
if ( dom(1) == dom(2) || dom(3) == dom(4) )
    error('ULTRASEM:ULTRASEMDOMAIN:rectangle', ...
        'Rectangle is degenerate.');
end
if ( nargin < 2 )
    m = 1;
end
if ( nargin < 3 )
    n = m;
end

% TODO: m and n are implemented the wrong way around below...
tmp = m;
m = n;
n = tmp;

% Initialize tree:
T.dom = dom;
T.idx = {};

if ( n == 1 && m == 1 )
    % Only a single patch. No subdivision or merge info required.
    T = ultraSEMDomain(T.dom, {[1 NaN]});
    return
end

x = linspace(dom(1), dom(2), m+1);
y = linspace(dom(3), dom(4), n+1);
l = 0;
for j = 1:n
    for k = 1:m
        l = l+1;
        T.dom(l,:) = [x(k), x(k+1), y(j), y(j+1)];
    end
end

if ( (log2(m) == round(log2(m))) && (log2(n) == round(log2(n))) )
    % Powers of two are easier to contend with

    l = 1;
    while ( m*n > 1 )
        tmp = reshape(1:m*n, m, n);
        if ( m >= n )
            % Merge up
            m = m/2;
        else
            % Merge right
            tmp = tmp.';
            n = n/2;
        end
        tmp = reshape(tmp, 2, m*n)';
        T.idx{l} = tmp; l = l + 1;
    end

else

    l = 0;
    while ( m*n > 1 )
        l = l+1;
        T.idx{l} = [];
        if ( n <= m )   % Merge right
            for k = 1:n
                idx = (k-1)*m + (1:m);
                if ( mod(m, 2) )
                    idx = [idx, NaN]; %#ok<AGROW>
                end
                idx = reshape(idx, 2, length(idx)/2)';
                T.idx{l} = [T.idx{l} ; idx];
            end
            m = ceil(m/2);
        else            % Merge up
            for k = 1:ceil(n/2)
                idx = 2*(k-1)*m + (1:2*m);
                idx(idx > m*n) = NaN;
                idx = reshape(idx, length(idx)/2, 2);
                idx(all(isnan(idx), 2),:) = [];
                T.idx{l} = [T.idx{l} ; idx];
            end
            n = ceil(n/2);
        end
    end

end

T = ultraSEMDomain(T.dom, T.idx);

end
