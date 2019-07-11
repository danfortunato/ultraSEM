function out = ultraSEMroot()
%ULTRASEMROOT   Root directory of ultraSEM installation.
%   S = ULTRASEMROOT() returns a string that is the name of the directory where
%   ultraSEM is installed.

out = fileparts(which('ultraSEMroot'));

end