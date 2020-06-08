function sol = solve(S, bc)
%SOLVE   Solve an ULTRASEM system.
%   SOL = SOLVE(S, BC) returns an ULTRASEM.SOL object SOL representing the
%   solution to the PDE, subject to the boundary conditions specified by
%   the function handle BC, contained in the ULTRASEM object S. If S has
%   not yet been built (see BUILD()) then SOLVE() will build it. If S has
%   not yet been initialized (see INITIALIZE()) then an error is thrown.
%
%   The full sequence for solving a problem using an ULTRASEM object S is:
%
%      initialize(S, OP, RHS)
%      build(S)
%      sol = S\bc % or sol = solve(S, bc)
%
%   See also BUILD, INITIALIZE.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

% Build the ULTRASEM if required:
if ( numel(S.patches) > 1 ), build(S); end

assert(isInitialized(S), ['The ultraSEM object `%s` has not been', ...
    'initialized.'], inputname(1));

% Build the boundary conditions:
coeffs = bc2coeffs(S, bc);

% Solve the patch object:
sol = solve(S.patches{1}, coeffs);

end
