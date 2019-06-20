function pass = test_mergeQuadsAndRects

% Merge a triangle with a square to make a house.

tol = 1e-8;
S = ultraSEM.rectangle([0 1 0 1], 2, 2);
T = ultraSEM.triangle([0 1 ; 1 1 ; .5 1+sqrt(3)/2]);

U = merge(S, T);
U = refine(U,1);

L = ultraSEM(U, {1,0,0}, -1, 21);
sol = L\0;

pass = abs(feval(sol, .6, .6) - 0.090868577361107) < tol;

end