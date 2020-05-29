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

assert(isInitialized(S), 'ultraSEM object `%s` has not been initialized.', ...
    inputname(1))

if ( isa(f, 'ultraSEM.Sol') )
    f = f.coeffs;
end

if ( iscell(f) )
    if ( isBuilt(S) )
        S.patches{1} = updateRHS(S.patches{1}, f);
    else
        for k = 1:numel(S.patches)
            S.patches{k} = updateRHS(S.patches{k}, f(k));
        end
    end
else
    for k = 1:numel(S.patches)
        S.patches{k} = updateRHS(S.patches{k}, f);
    end
end

end
