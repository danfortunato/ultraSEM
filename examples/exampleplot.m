function exampleplot(obj, k)

figure(1)
set(gcf, 'Position',  [0, 600, 1000, 400])

if ( isa(obj,'ultraSEM.Domain') )
    if ( nargin == 1 ), k = 1; end
    clf
    subplot(1,2,k)
    plot(obj)
    axis square tight
    grid off, box on
    subplot(122), axis off
    text(0, 0, 'Computing...', 'FontSize', 16, 'HorizontalAlignment', 'center')
    axis([-1 1 -1 1])
    shg
elseif ( isa(obj,'ultraSEM.Sol') )
    if ( nargin == 1 ), k = 2; end
    subplot(1,2,k), cla
    plot(obj)
    axis square tight
    grid off, box on
    pos = get(gca, 'Position');
    colorbar
    set(gca, 'Position', pos)
end

end