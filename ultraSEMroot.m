function out = ultraSEMroot()
%ULTRASEMROOT   Root directory of the ULTRASEM installation.
%   S = ULTRASEMROOT() returns a string that is the name of the directory
%   where ULTRASEM is installed.

out = fileparts(which('ultraSEMroot'));

end