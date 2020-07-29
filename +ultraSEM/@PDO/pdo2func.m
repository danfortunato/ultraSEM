function fh = pdo2func(L)
%PDO2CFUNC   Convert a PDO to a function handle.
%   FH = PDO2FUNC(L) creates a function handle representing the PDO L.
%
%   See also PDO2CELL.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

fh = [];
uniformOp = true;
for s = fieldnames(L).'
    s = s{1}; %#ok<FXSET>

    dx = length(strfind(s, 'x'));
    dy = length(strfind(s, 'y'));
    if ( dx == 0 && dy == 0 )
        term = 'u';
    else
        term = ['diff(u, [' num2str(dx) ' ' num2str(dy) '])'];
    end

    if ( isa(L.(s), 'function_handle') )
        uniformOp = false;
        coeff = regexprep(func2str(L.(s)), '@\(.+?\)', '');
        fh = [fh ' + (' coeff ').*' term]; %#ok<AGROW>
    elseif ( L.(s) ~= 0 )
        coeff = num2str(L.(s));
        fh = [fh ' + ' coeff '*' term]; %#ok<AGROW>
    end
end

% Strip off the leading +
if ( ~isempty(fh) )
    fh = fh(4:end);
else
    fh = '0*u';
end

if ( uniformOp )
    fh = ['@(u) ' fh];
else
    fh = ['@(x,y,u) ' fh];
end

fh = eval(fh);

end