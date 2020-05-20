function T = mtimes(T, c)
%*   Scale an ULTRASEM.SOL.
%   C*T will scale the ULTRASEM.SOL by C. C must be a scalar.

if ( ~isa(T, 'ultraSEM.Sol') )
    % Ensure T is the ULTRASEM.SOL:
    T = mtimes(c, T);
    return
elseif ( isa(c, 'ultraSEM.Sol' ) )
    % We can't multiply two ULTRASEM.SOLs:
    error('ULTRASEM:SOL:mtimes:twosols', ...
        'Cannot multiply (* or .*) two ultraSEM.Sols.\n')
elseif ( ~isnumeric(c) )
    error('ULTRASEM:SOL:mtimes:unknown', ...
        'Cannot multiply an object of type %s by an ultraSEM.Sol.', ...
        class(c));
elseif ( ~isscalar(c) )
    error('ULTRASEM:SOL:mtimes:scalar', ...
        'C must be a scalar.')
end

% Scale the solution:
T.u = cellfun(@(u) c*u, T.u, 'UniformOutput', false);

end
