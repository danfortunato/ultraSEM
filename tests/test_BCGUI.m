dom = ultraSEM.rectangle([-2 2 -1 1], 1, 2);
S = ultraSEM(dom, {1, 0, 0}, -1);
h = plot(S); shg

%%

bc = ultraSEM.BC.getBCsFromPlot(h);
plot(S\bc), shg

