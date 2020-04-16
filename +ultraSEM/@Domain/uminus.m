function dom = uminus(dom)
%-   Unary minus for an ULTRASEM.DOMAIN.
%   -DOM negates the ULTRASEM.DOMAIN DOM.
%
% See also UPLUS.

dom.domain = -dom.domain;

end
