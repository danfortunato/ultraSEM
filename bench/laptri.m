function t = laptri(r, p)

D = ultraSEM.Domain.triangle([0 0 ; 1 0 ; .5  sqrt(3)/2 ]);
D = refine(D, r);
op = {{1,0,1}, {0,0}, 0};

tic
for k = 1:1000  
    S = ultraSEM(D, op, -1, p);
    sol = S\0;
    t = toc;
    if ( t > 1 ), break, end
end
t = t/k;
    
end
