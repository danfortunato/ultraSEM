function pass = test_mergeL()

tol = 1e-10;

% Merges an L-shape and a square to form a square.

l = ultraSEM.alphabet('l');
s = ultraSEM.rectangle([0 1 0 1], 1);
D = l & s;
op = ultraSEM(D, {1, 0, 1}, -1, 21);
sol = op\0;

S = ultraSEM.rectangle([-1 1 -1 1], 2);
op = ultraSEM(S, {1, 0, 1}, -1, 21);
sol2 = op\0;

err = norm(sol - sol2, inf);
pass(1) = err < tol;

end
