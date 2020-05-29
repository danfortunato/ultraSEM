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
    [x, y] = getGrid(sol1);
    for j = 1:numel(sol1.coeffs)
        err = abs(feval(sol1, x{j}, y{j}) - feval(sol, x{j}, y{j}));
        pass(k) = pass(k) & max(max(err)) < tol;
    end
end

end
