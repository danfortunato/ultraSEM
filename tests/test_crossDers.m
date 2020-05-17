function pass = test_crossDers()

tol = 1e-10;

% Test that we solve for a constant function correctly when there are cross
% derivatives.

D = ultraSEM.rectangle([-1 1 -1 1]);
op = ultraSEM( D, {{1,1,1}, 0, 0}, 0, 51 );
sol = op \ 1;
% plot(sol)

err = norm(sol-1, inf);
pass = err < tol;

end
