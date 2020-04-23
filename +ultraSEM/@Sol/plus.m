function h = plus(f, g)
%+   Plus for ULTRASEM.SOL.
%   F + G adds the ULTRASEM.SOL objects F and G. F and G must have the same
%   domains and discretization sizes. F and G may also be scalars.
%
% See also MINUS.

if ( isnumeric( f ) )
    h = plus(g, f);
    return
elseif ( isnumeric( g ) )
    h = f;
    for k = 1:size(h.u,1)
        h.u{k}(1,1) = h.u{k}(1,1) + g;
    end
elseif ( isa(f, 'ultraSEM.Sol') && isa(g, 'ultraSEM.Sol') )
    h = f;
    % TODO: Assume on the same grid for now.
    h.u = cellfun(@plus, f.u , g.u, 'UniformOutput', false);
end

end