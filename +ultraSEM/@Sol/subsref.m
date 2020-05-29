function varargout = subsref(sol, index)
%SUBSREF   Subscripted reference for an ULTRASEM.SOL.
%   SOL(X, Y) returns the values of SOL evaluated at (X, Y). See FEVAL for
%   further details.
%
%   SOL.PROP returns the property PROP of SOL as defined by GET(SOL,
%   'PROP').
%
% See also FEVAL, GET.

idx = index(1).subs;
switch index(1).type

%% %%%%%%%%%%%%%%%%%%%%%%%%%% FEVAL / COMPOSE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case '()'

        % Where to evaluate:
        x = idx{1}; 
        y = idx{2}; 

        out = feval(sol, x, y);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GET %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case '.'

        % Call GET() for .PROP access.
        out = get(sol, idx);
        if ( numel(index) > 1 )
            % Recurse on SUBSREF():
            index(1) = [];
            out = subsref(out, index);
        end

    otherwise

        error('ULTRASEM:SOL:subsref:unexpectedType',...
            ['??? Unexpected index.type of ', index(1).type]);
end

% Convert to a cell:
varargout = {out};

end
