function dom = uminus(dom)
%-   Unary minus for an ULTRASEM.DOMAIN.
%   -DOM negates the ULTRASEM.DOMAIN DOM.
%
%   See also UPLUS.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

dom.domain = -dom.domain;

end
