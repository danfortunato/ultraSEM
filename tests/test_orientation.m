function pass = test_orientation()

tol = 1e-10;

sol = chebfun2(@(x,y) (1-x.^2).*(1-y.^2).*(x.^4.*y + y), [-1 1 -1 1]);
rhs = lap(sol);
bc = 0;
n = 10;
pdo = {{1, 0, 1}, {0, 0}, 0};
dom = [ -1 0  0 1 ;
         0 1  0 1 ;
        -1 0 -1 0 ;
         0 1 -1 0 ];

p = perms(1:4);
for k = 1:size(p,1)
    D = (ultraSEM.rectangle(dom(p(k,1),:))  & ...
         ultraSEM.rectangle(dom(p(k,2),:))) & ...
        (ultraSEM.rectangle(dom(p(k,3),:))  & ...
         ultraSEM.rectangle(dom(p(k,4),:)));
    op = ultraSEM(D, pdo, rhs, n);
    sol1 = op \ bc;
    pass(k) = 1;
    for j = 1:numel(sol1.u)
        xx = sol1.x{j}; yy = sol1.y{j};
        err = abs(feval(sol1, xx, yy) - feval(sol, xx, yy));
        pass(k) = pass(k) & max(max(err)) < tol;
    end
end

end
