function n = length(P)
%LENGTH   Number of patches in an ULTRASEM.PARENT.
%   LENGTH(P) returns the total number of patches in the ULTRASEM.PARENT
%   object P.

n = length(P.child1) + length(P.child2);

end
