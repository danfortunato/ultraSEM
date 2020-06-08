function [u, d] = solve(P, bc)
%SOLVE   Solve a leaf patch.
%   SOL = SOLVE(P, BC) returns an ULTRASEM.SOL object representing the PDE
%   solution on the leaf P with Dirichlet boundary data BC.
%
%   [U, D] = SOLVE(P, BC) returns instead a cell array containing the
%   solution coefficients U and a vector D containing the domain of the
%   patch.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

% Extract the domain from the patch:
d = P.domain;

% Evaluate the solution operator for the patch:
u = P.S * [bc ; 1]; % The 1 accounts for the particular part.
u = reshape(u, P.p*[1,1]);

% Return cell output for consistency with ULTRASEM.PARENT/SOLVE():
u = {u};

if ( nargout == 1 )
    % If a single output is requested, return an ULTRASEM.SOL object.
    u = ultraSEM.Sol(u, d);
end

end
