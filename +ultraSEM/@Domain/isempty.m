function out = isempty(T)
%ISEMPTY   Test for empty ULTRASEM.DOMAIN.
%   ISEMPTY(DOM) returns 1 if DOM is an empty ULTRASEM.DOMAIN and 0
%   otherwise.

out = isempty(T(1).domain);

end
