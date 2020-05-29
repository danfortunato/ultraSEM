function initialize(S, varargin)
%INITIALIZE   Initialize an ULTRASEM with a PDE and solve local subproblems.
%   INITIALIZE(S, OP, RHS) will initialize the ULTRASEM object S with the
%   PDE defined by OP and the righthand side RHS. The PDE will be solved on
%   each of the subpatches of S.domain.
%
%   INITIALIZE(S, OP) assumes the problem is homogeneous (i.e., RHS = 0).
%
%   INITIALIZE(S, OP, RHS, N) uses an N x N discretization on each patch.
%
%   INITIALIZE(..., PREF) uses the preferences specified in the
%   ULTRASEM.PREF object PREF. (See ULTRASEM.PREF for details on the
%   various preference options and their defaults.)
%
%   The full sequence for solving a problem using an ULTRASEM object S is:
%
%      initialize(S, OP, RHS)
%      build(S)   % (optional)
%      sol = S\bc % or sol = solve(S, bc)
%
% See also BUILD, SOLVE.

% Initialize all leaf patches:
D = S.domain;
numD = size(D.domain, 1);
if ( isa(D.domain, 'ultraSEM.Domain') && numD > 1)
    % Deal with recursively defined Domain:
    S_sub(numD,1) = ultraSEM;
    for k = 1:numD
        S_sub(k).domain = D.domain(k);
        initialize(S_sub(k), varargin{:});
    end
    S.patches = S_sub;
else
    S.patches = ultraSEM.Leaf.initialize(D.domain, varargin{:});
end

end
