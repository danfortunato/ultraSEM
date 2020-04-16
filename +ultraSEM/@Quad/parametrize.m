function obj = parametrize(obj)
%PARAMETRIZE   Compute the parametrization of an ULTRASEM.QUAD.

% Compute the change of variables:
M = [1 -1 -1  1;  % (-1,  1)
     1  1 -1 -1;  % ( 1, -1)
     1  1  1  1;  % ( 1,  1)
     1 -1  1 -1];
params = M \ obj.v;
obj.a1 = params(1,1); obj.a2 = params(1,2);
obj.b1 = params(2,1); obj.b2 = params(2,2);
obj.c1 = params(3,1); obj.c2 = params(3,2);
obj.d1 = params(4,1); obj.d2 = params(4,2);

end
