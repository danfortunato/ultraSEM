function out = ultraSEMroot()
%ULTRASEMROOT   Root directory of the ULTRASEM installation.
%   S = ULTRASEMROOT() returns a string that is the name of the directory
%   where ULTRASEM is installed.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

out = fileparts(which('ultraSEMroot'));

end