function t = gensquare(r, p)

D = ultraSEM.Domain.rectangle([-1 1 -1 1]);
D = refine(D, r);
op1 = @(x,y) cos(x+y);
op2 = @(x,y) 1./(1+exp(x-y));
op3 = 1;
op4 = @(x,y) sin(10*pi*x.*y);
op = {{op1,0,op2}, {op3,0}, op4};

tic
for k = 1:1000  
    S = ultraSEM(D, op, -1, p);
    sol = S\0;
    t = toc;
    if ( t > 1 ), break, end
end
t = t/k;
    
end
