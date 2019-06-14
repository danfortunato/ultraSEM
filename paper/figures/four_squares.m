%% Plot four squares in values

n = 7;
[x,y] = chebpts2(n);
bb = abs(x) == 1 | abs(y) == 1;
bb = bb(:);
ii = ~bb;

% Style
fc = 'b';
lc = 'r';
ms = 40;
lw = 4;

figure(1), clf
hold on

rectangle('Position',[0 0 2 2],'LineWidth',lw)
plot(x(ii)+1,y(ii)+1,'.','Color',fc,'MarkerSize',ms)
plot(x(bb)+1,y(bb)+1,'.','Color',lc,'MarkerSize',ms)

rectangle('Position',[-2 0 2 2],'LineWidth',lw)
plot(x(ii)-1,y(ii)+1,'.','Color',fc,'MarkerSize',ms)
plot(x(bb)-1,y(bb)+1,'.','Color',lc,'MarkerSize',ms)

rectangle('Position',[0 -2 2 2],'LineWidth',lw)
plot(x(ii)+1,y(ii)-1,'.','Color',fc,'MarkerSize',ms)
plot(x(bb)+1,y(bb)-1,'.','Color',lc,'MarkerSize',ms)

rectangle('Position',[-2 -2 2 2],'LineWidth',lw)
plot(x(ii)-1,y(ii)-1,'.','Color',fc,'MarkerSize',ms)
plot(x(bb)-1,y(bb)-1,'.','Color',lc,'MarkerSize',ms)

hold off
axis([-2.1 2.1 -2.1 2.1])
axis square off

matlab2tikz('values.tex')

%% Plot four squares in coefficients

% Style
fc = 'b';
lc = 'r';
lw = 6;
pad = 0.2;

figure(1), clf
axis tight square off
hold on

rectangle('Position',[pad pad 2 2],'FaceColor',fc)
rectangle('Position',[-2-pad pad 2 2],'FaceColor',fc)
rectangle('Position',[pad -2-pad 2 2],'FaceColor',fc)
rectangle('Position',[-2-pad -2-pad 2 2],'FaceColor',fc)

line([0 0],[ 0 2]+pad,'Color',lc,'LineWidth',lw)
line([0 0],[-2 0]-pad,'Color',lc,'LineWidth',lw)
line([-2 0]-pad,[0 0],'Color',lc,'LineWidth',lw)
line([ 0 2]+pad,[0 0],'Color',lc,'LineWidth',lw)

line([2 2] + 2*pad,[ 0 2]+pad,'Color',lc,'LineWidth',lw)
line([2 2] + 2*pad,[-2 0]-pad,'Color',lc,'LineWidth',lw)
line([-2 -2]-2*pad,[ 0 2]+pad,'Color',lc,'LineWidth',lw)
line([-2 -2]-2*pad,[-2 0]-pad,'Color',lc,'LineWidth',lw)

line([-2 0]-pad,[2 2]+2*pad,'Color',lc,'LineWidth',lw)
line([ 0 2]+pad,[2 2]+2*pad,'Color',lc,'LineWidth',lw)
line([-2 0]-pad,[-2 -2]-2*pad,'Color',lc,'LineWidth',lw)
line([ 0 2]+pad,[-2 -2]-2*pad,'Color',lc,'LineWidth',lw)

hold off
axis([-2.5 2.5 -2.5 2.5])
axis square off

matlab2tikz('coeffs.tex')
