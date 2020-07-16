function varargout = subsref(S, index)
%SUBSREF   Subscripted reference of an ULTRASEM object.
%   S(U) returns the forward application of the ULTRASEM S to the
%   ULTRASEM.SOL U. This is equivalent to S*U.
%
%   See also MTIMES.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

idx = index(1).subs;

switch ( index(1).type )

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% APPLY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case '()'

        % Forward application of an ULTRASEM to an ULTRASEM.SOL.
        if ( length(idx) == 1 && isa(idx{1}, 'ultraSEM.Sol') )
            % Call MTIMES:
            u = idx{1};
            out = S*u;
        else
            error('ULTRASEM:ULTRASEM:mtimes:badInput', ...
                'Can only forward apply an ULTRASEM to a single ULTRASEM.SOL.');
        end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GET %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case '.'

        % Call GET() for .PROP access.
        out = get(S, idx);
        if ( numel(index) > 1 )
            % Recurse on SUBSREF():
            index(1) = [];
            out = subsref(out, index);
        end

    otherwise

        error('ULTRASEM:ULTRASEM:subsref:unexpectedType', ...
            ['??? Unexpected index.type of ', index(1).type]);

end

% Convert to a cell:
varargout = {out};

end
