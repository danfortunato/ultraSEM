function [op, rhs] = transformPDO(dom, op, rhs)
%TRANSFORMPDO   Convert PDO on an ULTRASEM.RECT to equivalent PDO on [-1,1]^2.

if ( nargin > 2 && isa(rhs, 'function_handle') )
    rhs = @(r,s) rhs(dom.x(r,s), dom.y(r,s));
end

end
