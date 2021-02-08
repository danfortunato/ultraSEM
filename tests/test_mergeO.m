function pass = test_mergeO()

tol = 1e-10;

% Merges an O-shape and a square to form a larger square.

o = ultraSEM.alphabet('o');
s = ultraSEM.rectangle([1 2 1 2]);
D = (o & s) - (1+1i);
op = ultraSEM(D, {1, 0, 1}, -1, 21);
sol = op\0;

S = ultraSEM.rectangle([-1 2 -1 2], 3, 3);
op = ultraSEM(S, {1, 0, 1}, -1, 21);
sol2 = op\0;

err = abs(feval(sol,.5,.5) - feval(sol2,.5,.5));
pass(1) = err < tol;

end
