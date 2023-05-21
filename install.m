function success = install(path)
%INSTALL installs MatCUTEst at fullfile(path, 'matcutest'), using the compiled CUTEst in the
% directory containing this M file or at zipurl.

if ispc || ismac
    error('MatCUTEst:InvalidOS', 'This package supports only GNU/Linux.');
end

% 7z is needed. For Ubuntu, try `type 7z || sudo apt install p7zip-full`.
[status, zpath] = system('command -v 7z');
if status ~= 0 || isempty(zpath)
    error('MatCUTEst:SevenzNotFound', 'Command ''7z'' not found. Make sure that p7zip-full is installed.');
end

% pkg is the package containing the compiled version of MatCUTEst.
pkg = 'full.matcutest.7z.001';
gitaccount = 'equipez';
gitrepo = 'matcutest_compiled';
zipurl = ['https://github.com/', gitaccount, '/', gitrepo, '/archive/refs/heads/main.zip'];

% mdir is the directory containing this script.
mdir = fileparts(mfilename('fullpath'));
% cpwd is the current directory.
cpwd = pwd();

if nargin == 0
    % Decide `path` if it is not provided.
    homedir = getenv('HOME');
    if exist(fullfile(homedir, 'local'), 'dir')
        path = fullfile(homedir, 'local');
    elseif exist(fullfile(homedir, '.local'), 'dir')
        path = fullfile(homedir, '.local');
    else
        path = fileparts(mdir);  % The parent directory of mdir
    end
else
    if ~exist(path, 'dir')
        mkdir(path);
    end
end

% matcutest is the root directory of MatCUTEst.
matcutest = fullfile(path, 'matcutest');
if exist(matcutest, 'dir') || exist(matcutest, 'file')
    fprintf('\nThe package path\n\n    %s\n\nalready exists. Remove it to (re-)install the package.\n\n',  matcutest);
    success = false;
    return
end
fprintf('\nMatCUTEst will be installed at\n\n    %s\n', matcutest);

if exist(fullfile(mdir, pkg), 'file')
    pkg = fullfile(mdir, pkg);
else
    [~, name, ext] = fileparts(zipurl);
    zipname = [gitrepo, '-', name, ext];
    zipname = fullfile(tempdir, zipname);

    fprintf('\nDownloading the compiled package from\n\n    %s\n', zipurl);
    websave(zipname, zipurl);
    unzip(zipname, tempdir);
    zipdir = fullfile(tempdir, [gitrepo, '-', name]);
    pkg = fullfile(zipdir, pkg);
    fprintf('\nDone.\n');
end

fprintf('\nInstalling the compiled package\n\n    %s\n', pkg);

exception = [];
try
    cd(matcutest);
    dir
    cd(path);
    % Run `7z x pkg`. Use `evalc` to make is quiet.
    evalc('system([''7z x '', pkg])');
    cd(fullfile(matcutest, 'mtools'));
    setup();
catch exception
    % Do nothing for the moment.
end
cd(cpwd);

if isempty(exception)
    fprintf('MatCUTEst is successfully installed at\n\n    %s\n', matcutest);
    fprintf('\nTry "help matcutest" for more information.\n\n');
    success = true;
else
    rethrow(exception);
end

return
