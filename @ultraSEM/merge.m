function h = merge(f, g)
%MERGE   Merge two ULTRASEM objects.
%   H = MERGE(F, G) will merge the two ULTRASEM objects F and G. If F and G
%   have already been initialized and built, then H will be also.
%
% See also INITIALIZE, BUILD.

% Merge the domains:
h = ultraSEM();
h.domain = merge(f.domain, g.domain);

% Merge the patches:
if ( numel(f.patches) == 1 && numel(g.patches) == 1)
    % f and g have already been built.
    h.patches{1} = merge(f.patches{1}, g.patches{1});
else
    h.patches = [f.patches ; g.patches];
    h.initialize;
end

end
