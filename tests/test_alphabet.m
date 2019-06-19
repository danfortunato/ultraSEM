function pass = test_alphabet()

tol = 1e-13;
n = 20;

% Check solver on connected domain with the real part of holomorphic function:
S = ultraSEMDomain( 'S' ); 
f = @(x,y) real(exp(-(x+1i*y))); 
op = ultraSEM(S, {{1,0,1}, 0, 0}, 0, n); 
sol = op \ f;
[x, y] = getGrid(sol, 1);
pass(1) = ( norm(feval(sol,x,y) - f(x,y)) ) < tol; 

% Check solver on nonconnected domain with the real part of holomorphic function:
S = ultraSEMDomain( 'U' ); 
f = @(x,y) real(exp(-(x+1i*y))); 
op = ultraSEM(S, {{1,0,1}, 0, 0}, 0, n); 
sol = op \ f;
[x, y] = getGrid(sol, 1);
pass(2) = ( norm(feval(sol,x,y) - f(x,y)) ) < tol; 

end
