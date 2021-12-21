function pass = test_integral2()

tol = 1e-15;

dom = ultraSEM.rectangle([0 1 0 1]);
area = 1;
u = ultraSEM.Sol(@(x,y) 1+0*x, 5, dom);
pass(1) = ( norm(sum2(u) - area) ) < tol;

dom = ultraSEM.polygon(4);
area = 2;
u = ultraSEM.Sol(@(x,y) 1+0*x, 5, dom);
pass(2) = ( norm(sum2(u) - area) ) < tol;

dom = ultraSEM.polygon(5);
a = 2*sin(4*pi/5);
area = a^2 * sqrt(5*(5+2*sqrt(5)))/4;
u = ultraSEM.Sol(@(x,y) 1+0*x, 5, dom);
pass(3) = ( norm(sum2(u) -  area) ) < tol;

dom = ultraSEM.polygon(6);
area = 3*sqrt(3)/2;
u = ultraSEM.Sol(@(x,y) 1+0*x, 5, dom);
pass(4) = ( norm(sum2(u) - area) ) < tol;

end
