function pass = test_square

pass = true;
return

tol = 1e-6;

R = ultraSEM.rectangle;
T = ultraSEM.triangle;
T2 = ultraSEM.triangle([1 1 ; 2 1 ; 1.5 1-sqrt(3)/2]);

pass = true(1,3);

D = .5*R & (.5*R+2);
try
    op = ultraSEM(D, {1,0,0}, -1);
    u = op\0;
catch
    pass(1) = false;
end

D = T & T2;
try
    op = ultraSEM(D, {1,0,0}, -1);
    u = op\0;
catch
    pass(2) = false;
end

D = T & (1/2.5*R-.5+.4i);
try
    op = ultraSEM(D, {1,0,0}, -1);
    u = op\0;
catch
    pass(3) = false;
end

end

