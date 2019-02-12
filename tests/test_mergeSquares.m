function pass = test_mergeSquares()

tol = 1e-10;

% Test that we can merge [A, B, C] -> [A, B, A] -> [A A A]
%  (i.e., that we can merge domains which don't intersect).

s1 = ultraSEM.rectangle([0 1 0 1]);
s2 = ultraSEM.rectangle([1 2 0 1]);
s3 = ultraSEM.rectangle([2 3 0 1]);
S = (s1 & s2) & s3;
op = ultraSEM(S, {1, 0, 1}, -1, 21);
sol = op\0;

S2 = ultraSEM.rectangle([0 3 0 1], 1, 3);
op = ultraSEM(S2, {1, 0, 1}, -1, 21);
sol2 = op\0;

err(1) = abs(feval(sol,1.5,.5) - feval(sol2,1.5,.5));

%%

S3 = (s1 & s3) & s2;
op = ultraSEM(S3, {1, 0, 1}, -1, 21);
sol3 = op\0;

err(2) = abs(feval(sol,1.5,.5) - feval(sol3,1.5,.5));

%%

pass = err < tol;

end
