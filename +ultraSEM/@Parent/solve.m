function [u, d] = solve(P, bc)
%SOLVE   Solve a parent patch.
%   SOL = SOLVE(P, BC) returns an ULTRASEM.SOL object representing the PDE
%   solution on the parent P with Dirichlet boundary data given by BC.
%
%   [U, D] = SOLVE(P, BC) returns instead cell arrays containing the
%   solution coefficients U and a vector D containing the domains of the
%   subpatches.

% Evaluate the solution operator for the parent:
u = P.S * [bc ; 1]; % The 1 accounts for the particular part.

% Construct boundary conditions for children and recurse.

% Construct boundary indices:
i1 = 1:numel(P.idx1{1});
i2 = 1:numel(P.idx2{1});
if ( ~isempty(i1) )
    i2 = i2 + i1(end);
end
% CAT() is 10x faster than CELL2MAT().
idx1 = cat(1, P.idx1{:}); % idx1 = cell2mat(p.idx1.');
idx2 = cat(1, P.idx2{:}); % idx2 = cell2mat(p.idx2.');

% Assemble boundary conditions for child patches:
% - Apply the global-to-local map (i.e., the transpose of the
%   local-to-global map) to obtain the local boundary conditions.
ubc1 = ones(size(P.child1.S, 2)-1, 1);
ubc1(idx1) = [bc(i1) ; P.l2g1.'*u];
ubc2 = ones(size(P.child2.S, 2)-1, 1);
ubc2(idx2) = [bc(i2) ; P.l2g2.'*u];

% Solve for the child patches:
[u1, d1] = solve(P.child1, ubc1);
[u2, d2] = solve(P.child2, ubc2);

% Concatenate for output:
u = [u1 ; u2]; 
d = [d1 ; d2];

if ( nargout == 1 )
    % If a single output is requested, return an ULTRASEM.SOL object.
    u = ultraSEM.Sol(u, d);
end

end
