% N.B.:
% 1. This file should be placed at the userpath, which is by default
% fullfile(matlabroot,'toolbox', 'local')
% To see whether the userpath is set to the above, run
% userpath
% 2. If the userpath is set as the above, then do
% ln -s $HOME/Documents/iUbuntu/MATLABROOT%toolbox%local%startup.m $MATLABROOT/toolbox/local/startup.m
% where $MATLABROOT should be the value of `matlabroot`.
% 3. It may happen (strangely) that the userpath is unset (empty). In that case, run
% userpath(fullfile(matlabroot,'toolbox', 'local')); rehash toolboxcache
% where "rehash toolboxcache" updates the toolbox cache. Without updating the cache, MATLAB may not
% see the newly added file under toolbox/local. See
% https://www.mathworks.com/help/matlab/matlab_env/toolbox-path-caching-in-the-matlab-program.html


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Default setting from MATLAB R2022b
%STARTUPSAV   Startup file
%   Change the name of this file to STARTUP.M. The file
%   is executed when MATLAB starts up, if it exists
%   anywhere on the path.  In this example, the
%   MAT-file generated during quitting using FINISHSAV
%   is loaded into MATLAB during startup.

%   Copyright 1984-2000 The MathWorks, Inc.

%load matlab.mat  % Does not exist???
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


homedir = getenv('HOME');
logdir = fullfile(homedir, 'local', 'log', 'matlab');
if ~exist(logdir, 'dir')
    mkdir(logdir);
end

% Set the path for saving crash dumps and logs.
%setenv('MATLAB_LOG_DIR', logdir);  % Seems not working!?
% The following suggested by GitHub copilot does not work either.
%com.mathworks.services.Prefs.setBooleanPref('CrashDumpsSaveToHomeDir', false);
%com.mathworks.services.Prefs.setStringPref('CrashDumpsSaveDir', logdir);

% Turn on diary
diaryname = sprintf('matlab_%s.log', datestr(now, 'yymmdd_HHMMSS'));
diaryfile = fullfile(logdir, diaryname);
if exist(diaryfile, 'file')
    delete(diaryfile);
end
diary(diaryfile);
fprintf('The diary of this session is saved at\n\n\t%s\n\n', diaryfile);
% diary on;  % No need, as diary(diaryfile) turns on diary.

% Prevent the MATLAB Editor/Debugger from opening a MATLAB file when encountering a breakpoint
com.mathworks.services.Prefs.setBooleanPref('EditorGraphicalDebugging',false);

% Set default Figure Position.
sizeofscreen = get(0, 'ScreenSize');
set(0, 'DefaultFigurePosition', [sizeofscreen(3)*0.7, sizeofscreen(4)-0.93*sizeofscreen(3)*0.3, sizeofscreen(3)*0.3, 0.81*sizeofscreen(3)*0.3]);
