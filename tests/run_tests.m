function run_tests

testFiles = dir('test_*.m');
testFiles = {testFiles.name};

% For making the output string align nicely:
maxLength = max(cellfun(@length, testFiles));

% Allocate pass and timing variables:
numFiles  = numel(testFiles);
durations = zeros(numFiles, 1);
errorMessages = {'FAILED', 'CRASHED'};

testDir = '';

% Attempt to run all of the tests:
try % Note, we try-catch as we've CD'd and really don't want to end up elsewhere
    
    % Loop over the test files:
    for k = 1:numFiles
        % Next file to test: (.m extension is removed).
        testFile = testFiles{k}(1:end-2);
       

        printTestInfo(testDir, testFile, k, maxLength);
        durations(k) = runTest(testFile);
        if ( durations(k) > 0 )
            % Success message.
            message = sprintf('passed in %.4fs', durations(k));
        else
            % Error flags are negative indices.
            message = errorMessages{-durations(k)};
        end
        fprintf([message '\n']);

    end
    
catch ME
    
    % We failed. Return to the starting directory and rethrow the error:
    rethrow(ME)
    
end

% Yay!
if ( all(durations > 0) )
    fprintf('All %s tests passed in %4.4fs.\n', testDir, sum(durations));
elseif ( ~quietMode )
    % Note. We don't show this in quiet mode as it's already clear.
    fprintf('%d failed test(s) in %s directory.\n', sum(durations < 0), testDir);
end

end


function printTestInfo(testDir, testFile, k, maxLength)
%PRINTTESTINFO   Pretty print test info for a given test file.
%   PRINTTESTINFO(TESTDIR, TESTFILE, K) will print the following to screen:
%     Test #k: TESTDIR/TESTFILE.m ...
%   Note that it will not linebreak after this.
%
%   PRINTTESTINFO(TESTDIR, TESTFILE, K, MAXLENGTH) is similar, but will print
%   extra whitespace after the ellipsis so that the following text for each
%   TESTFILE in TESTDIR is aligned.
%
%   Furthermore, PRINTTESTINFO tests to see if the Matlab desktop is running and
%   will support HTML. If it is, the displayed the displayed TESTDIR/TESTFILE.m
%   is hyperlinked to open TESTFILE.m in the Matlab editor.

link = printTestFilename(testDir, testFile);           % Either html or plain.
ws = repmat(' ', 1, maxLength - length(testFile) - 1); % Whitespace.
fprintf('  Test #%3.3d: %s ...%s', k, link, ws);       % Print to screen.

end

function link = printTestFilename(testDir, testFile, loc)
%PRINTTESTFILENAME   Pretty print a filename.
%   STR = PRINTTESTFILENAME(TESTDIR, TESTFILE, LOC) tests to see if the Matlab
%   desktop is running and will support HTML. If it is, the displayed the
%   displayed TESTDIR/TESTFILE.m is hyperlinked to open LOC/TESTFILE.m in the
%   Matlab editor. If LOC is not passed, PRINTTESTFILENAME will attempt to find
%   the location of TESTFILE.M via fileparts(which(TESTFILE)).

if ( usejava('jvm') && usejava('desktop') )
    % Use HTML links if Java is enabled. 
    
    if ( nargin < 3 )
        loc = fileparts(which(testFile));
    end
    url = fullfile(loc, [testFile '.m']);
    
    link = ['<a href = "matlab: edit ''', url, '''">' testDir '/', ...
        testFile '.m</a>'];  
    
else
    
    % Otherwise, use plaintext.
    link = fullfile(testDir, [testFile '.m']);
    
end

end


function duration = runTest(testFile)
%RUNTEST Runs the given test file.
%   DURATION = RUNTEST(TESTFILE) executes the file TESTFILE.
%   TESTFILE should return a vector or matrix of logical values. 
%
%   If each of these are logical true (or 1) then the test passes and RUNTEST
%   returns DURATION as a double > 0 corresponding to the time it took the test
%   to execute in seconds.
%
%   If any of the entries returned by TESTFILE are logcial false (or 0), then
%   DURATION = -1.
%
%   If executing TESTFILE crashes, this is caught in a try-catch statement, and
%   RUNTTEST returns DURATION = -2.

% Store current default preference states:
prefState1 = chebfunpref();
prefState2 = cheboppref();
% Close any open windows:
close all

% Attempt to run the test:
try
    tstart = tic();
    pass = feval(testFile);
    duration = toc(tstart);

    pass = all(pass(:));
    if ( ~pass )
        % The test failed, so return FAILED flag.
        duration = -1;
    end

catch ME %#ok<NASGU>
    % We crashed. This is bad. Return CRASHED flag.
    duration = -2;
    
    % But we _don't_ want to throw an error.
    %rethrow(ME)
end

% Close any windows the test may have left open:
close all
% Ensure global preferences aren't modified by tests.
chebfunpref.setDefaults(prefState1);
cheboppref.setDefaults(prefState2);

end