function pass = test_mergeQuadsAndRects

% Merge a triangle with a square to make a house.

tol = 1e-6;
S = ultraSEM.rectangle([0 1 0 1], 2, 2);
T = ultraSEM.triangle([0 1 ; 1 1 ; .5 1+sqrt(3)/2]);

U = merge(S, T);
U = refine(U,1);

L = ultraSEM(U, {1,0,0}, -1, 15);
sol = L\0;

pass = abs(feval(sol, .6, .6) - 0.090868577361107) < tol;

% Nick's bug: Merging a Quad and a Rect was wrong if the Quad was
% initialized before the Rect. (Note that merging an ultraSEM.quad with an
% ultraSEM.rectangle was fine---it was specifically the Quad and Rect
% objects.)
D = ultraSEM.quad([-1 0; -1.5 1; 0 1; 0 0]) & ultraSEM.rectangle([0 1 0 1]);
D.domain = cat(1, D.domain.domain); % Convert to Quad and Rect
S = ultraSEM(D, {1, 0, @(x,y) 10*(x.^2+y.^2)}, -1);
sol = S \ 0;

% Now force both to be Quad
D.domain(2) = ultraSEM.Quad(quadVertices(D.domain(2)));
S = ultraSEM(D, {1, 0, @(x,y) 10*(x.^2+y.^2)}, -1);
sol_true = S \ 0;

pass(2) = norm(sol - sol_true, inf) < tol;

end