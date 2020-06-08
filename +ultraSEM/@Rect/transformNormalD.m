function normal_d = transformNormalD(T, p)
%TRANSFORMNORMALD   Normal derivative operator for rectangular domains.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

persistent bcrows_d % Store for efficiency.
if ( size(bcrows_d, 2) ~= p )
    % Construct normal derivatives along the four edges:
    I = speye(p);
    lbc_d = kron( (-1).^(0:p-1).*(0:p-1).^2, I );
    rbc_d = kron( ones(1,p).*(0:p-1).^2, I );
    dbc_d = kron( I, (-1).^(0:p-1).*(0:p-1).^2 );
    ubc_d = kron( I, ones(1,p).*(0:p-1).^2 );
    bcrows_d = [ lbc_d ; rbc_d ; dbc_d ; ubc_d ];
end

rect = rectVertices(T);
domx = rect(1:2);    domy = rect(3:4); 
sclx = 2/diff(domx); scly = 2/diff(domy);
normal_d = bcrows_d;
normal_d(1:2*p,:) = sclx*normal_d(1:2*p,:);
normal_d(2*p+1:end,:) = scly*normal_d(2*p+1:end,:);

end
