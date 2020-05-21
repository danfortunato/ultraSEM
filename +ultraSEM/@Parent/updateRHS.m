function P = updateRHS(P, rhs)
%UPDATERHS   Update RHS of an ULTRASEM.PARENT object.
%   P = UPDATERHS(P, RHS) replaces the existing RHS of an initialized
%   ULTRASEM.PARENT object P with that given in RHS, which may be a
%   constant, a function handle, an ULTRASEM.SOL object defined on P, or a
%   cell array containing bivariate Chebyshev coefficients for each patch.

if ( isa(rhs, 'ultraSEM.Sol') )
    rhs = rhs.u;
end

if ( iscell(rhs) )
    % We are updating the RHS from a cell array of coefficients.
    % Get the number of patches in each child:
    n1 = length(P.child1);
    n2 = length(P.child2);
    % Update RHS of children:
    a = updateRHS(P.child1, rhs(1:n1));
    b = updateRHS(P.child2, rhs(n1+1:n1+n2));
else
    % Update RHS of children:
    a = updateRHS(P.child1, rhs);
    b = updateRHS(P.child2, rhs);
end

i1 = P.idx1{1};
s1 = P.idx1{2};
i2 = P.idx2{1};
s2 = P.idx2{2};
l2g1 = P.l2g1;
l2g2 = P.l2g2;

% Extract D2N maps:
D2Na = a.D2N; D2Nb = b.D2N;
% and discard from children
% a.D2N = []; b.D2N = [];

% Compute new solution operator:
S = solve(P.dA, l2g1*D2Na(s1,end) + l2g2*D2Nb(s2,end), false);
%              |------------------ rhs --------------|

% Compute new D2N maps:
%      |--- rhs ----|
D2N = [ D2Na(i1,end) ;
        D2Nb(i2,end) ] ...
    + [ D2Na(i1,s1)*l2g1.' ; D2Nb(i2,s2)*l2g2.' ] * S;

P.S(:,end) = S;
P.D2N(:,end) = D2N;

P.child1 = a;
P.child2 = b;

end
