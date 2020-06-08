function pref = mergePrefStructs(pref1, pref2)
%MERGEPREFSTRUCTS   Merge preference structures.
%   P = MERGEPREFSTRUCTS(P, Q), where P and Q are MATLAB structures,
%   "merges" Q into P by replacing the contents of fields in P with those
%   of identically-named fields in Q. If Q has a field whose name does not
%   match any of those in P, it is added to P.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

pref = pref1;
for field = fieldnames(pref2).'
    pref.(field{1}) = pref2.(field{1});
end

end
