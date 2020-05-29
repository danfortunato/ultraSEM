function S = updateRHS(S, f)
%UPDATERHS   Update the RHS function in a ULTRASEM object.
%   UPDATERHS(S, F) or S.RHS = F updates the RHS function F in the ULTRASEM
%   object S. The function F may be a constant, a function handle, an
%   ULTRASEM.SOL object defined on S, or a cell array containing bivariate
%   Chebyshev coefficients for each patch. This is useful as only a small
%   fraction of the initialization phase needs to be repeated.
%
% Example:
%   S = ultraSEM.alphabet('S');
%   T = refine(S, 2);
%   op = {1, 0, @(x,y) y};
%   tic
%       S = ultraSEM(T, op, -1);
%       S.build;
%       sol1 = S\0;
%   t1 = toc;               % t1 = 0.866088
%   tic
%       S.rhs = @(x,y) sin(pi*x.*y);
%       sol2 = S\0;
%   t2 = toc;               % t2 = 0.277288
%
% See also BUILD.

assert(isInitialized(S), 'ultraSEM object has not been initialized.')

% Extract coeffs if f is a ultraSEM.Sol:
if ( isa(f, 'ultraSEM.Sol') )
    f = f.coeffs;
end
% Duplicate f if it is a scalar or function handle:
if ( ~iscell(f) )
    f = repmat({f}, numel(S.patches), 1);
end
% Update RHS of each patch:
for k = 1:numel(S.patches)
    S.patches{k} = updateRHS(S.patches{k}, f{k});
end

end
