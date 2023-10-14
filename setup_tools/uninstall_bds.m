function uninstall_bds(path_string_stamp)
%This file is cited from https://github.com/libprima/prima/blob/main/matlab/setup_tools/uninstall_prima.m,
%which is written by Zaikun Zhang.
%UNINSTALL_BDS uninstalls BDS.

fprintf('\nUninstalling BDS (if it is installed) ... ');

% The full path of several directories.
mfiledir = fileparts(mfilename('fullpath'));  % The directory where this .m file resides
matd = fileparts(mfiledir); % Matlab directory
src = fullfile(matd, 'src'); % Directory containing the source code
examples = fullfile(matd, 'examples'); % Directory containing some examples
tests = fullfile(matd, 'tests'); % Directory containing some tests
competitors = fullfile(tests, 'competitors'); % Directory containing some competitors
tools = fullfile(matd, 'setup_tools'); % Directory containing some tools
path_string = {src, examples, tests, competitors}; % The paths to be removed

% Try removing the paths possibly added by PRIMA
orig_warning_state = warning;
warning('off', 'MATLAB:rmpath:DirNotFound'); % Maybe the paths were not added. We do not want to see this warning.
warning('off', 'MATLAB:SavePath:PathNotSaved'); % Maybe we do not have the permission to save path.
rmpath(src, examples, tests, competitors, tools);
savepath;
warning(orig_warning_state); % Restore the behavior of displaying warnings

% Removing the line possibly added to the user startup script
user_startup = fullfile(userpath,'startup.m');
if exist(user_startup, 'file')
    for i = 1:length(path_string) 
        add_path_string = sprintf('addpath(''%s'');', path_string{i});
        full_add_path_string = sprintf('%s\t%s %s', add_path_string, '%', path_string_stamp);
        try
            del_str_ln(user_startup, full_add_path_string);
        catch
            % Do nothing.
        end
    end
end

callstack = dbstack('-completenames');
root_dir = fileparts(callstack(2).file);  % Root directory of the package
fprintf('Done.\nYou may now remove\n\n    %s\n\nif it contains nothing you want to keep.\n\n', root_dir);

return
