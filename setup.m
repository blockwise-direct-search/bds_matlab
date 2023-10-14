function setup(varargin)
%This file is cited from https://github.com/libprima/prima/blob/main/setup.m, which is written by Zaikun Zhang.
%SETUP sets the package up for MATLAB.
%
%   This script can be called in the following ways.
%
%   setup % Add the paths needed to use the package
%   setup solvers % Add the paths needed to use the solvers
%   setup uninstall  % Uninstall the package
%
%   REMARKS:
%
%   1. To run this script, you need to have write access to the directory that
%   contains this script and its subdirectories.
%
%   ***********************************************************************
%   Authors:    Haitian LI (hai-tian.li@connect.polyu.hk)
%               and Zaikun ZHANG (zaikun.zhang@polyu.edu.hk)
%               Department of Applied Mathematics,
%               The Hong Kong Polytechnic University.
%   ***********************************************************************

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Attribute: public (can be called directly by users)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% setup starts

% path_stringName of the package. It will be used as a stamp to be included in the path_string. Needed only
% if `savepath` fails.
package_name = 'bds';

% Check the version of MATLAB.
if verLessThan('matlab', '9.2')   % MATLAB R2017a = MATLAB 9.2
    fprintf('\nSorry, this package does not support MATLAB R2017b or earlier releases.\n\n');
    return
end

% The full paths to several directories needed for the setup.
setup_dir = fileparts(mfilename('fullpath')); % The directory containing this setup script.
setup_tools_dir = fullfile(setup_dir, 'setup_tools'); % Directory containing some tools for setting up.
src_dir = fullfile(setup_dir, 'src'); % Directory containing the source code of the package.
examples_dir = fullfile(setup_dir, 'examples'); % Directory containing some examples.
tests_dir = fullfile(setup_dir, 'tests'); % Directory containing some tests
tests_competitors_dir = fullfile(tests_dir, 'competitors'); % Directory containing some tests for competitors.


% We need write access to `setup_dir` (and its subdirectories). Return if we do not have it.
% N.B.: This checking is NOT perfect because of the following --- but it is better than nothing.
% 1. `fileattrib` may not reflect the attributes correctly, particularly on Windows. See
% https://www.mathworks.com/matlabcentral/answers/296657-how-can-i-check-if-i-have-read-or-write-access-to-a-directory
% 2. Even if we have write access to `setup_dir`, we may not have the same access to its subdirectories.
[~, attribute] = fileattrib(setup_dir);
if ~attribute.UserWrite
    fprintf('\nSorry, we cannot continue because we do not have write access to\n\n%s\n\n', setup_dir);
    return
end

% `tools` contains some functions needed in the sequel.
addpath(setup_tools_dir);

% Parse the input.
[action, wrong_input] = parse_input(varargin);

% Exit if wrong input detected. Error messages have been printed during the parsing.
if wrong_input
    rmpath(setup_tools_dir);
    error('bds:InvalidInput', 'setup: The input is invalid.');
end

% Uninstall the package if requested.
if strcmp(action, 'uninstall')
    uninstall_bds(package_name);
    return
end


%%%%%%%%%%%%%%% If we arrive here, then the user requests us to compile the solvers. %%%%%%%%%%%%%%%

% Add the related paths to the MATLAB path, and then try saving the path.
if strcmp(action, 'bds')
    path_saved = add_save_path({src_dir}, package_name);
else 
    path_saved = add_save_path({src_dir, examples_dir, tests_dir, tests_competitors_dir}, package_name);
end

fprintf('\nThe package is ready to use.\n');
fprintf('\nYou may now try ''help bds'' for information on the usage of the package.\n');

if ~strcmp(action, 'bds')
    fprintf('\nYou may also run ''testbds'' to test the package on a few examples.\n');
end

rmpath(setup_tools_dir);

if ~path_saved  % `add_save_path` failed to save the path.
    add_path_string = sprintf('addpath(''%s'');', src_dir);
    fprintf('\n***** To use the package in other MATLAB sessions, do ONE of the following. *****\n');
    fprintf('\n- Append the following line to your startup script');
    fprintf('\n  (see https://www.mathworks.com/help/matlab/ref/startup.html for information):\n');
    fprintf('\n    %s\n', add_path_string);
    fprintf('\n- OR come to the current directory and run ''setup'' when you need the package.\n');
end

fprintf('\n');

% setup ends
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function path_saved = add_save_path(path_string, path_string_stamp)
%ADD_SAVE_PATH adds the path indicated by PATH_STRING to the MATLAB path and then tries saving path.
% PATH_STRING_STAMP is a stamp used when writing PATH_STRING to the user's startup.m file, which is
% needed only if `savepath` fails.
% N.B.: Why not putting this function as an individual file in the `tools` directory? Because we
% need it after `tools` is removed from the path.

if nargin < 2
    path_string_stamp = sprintf('Added by %s', mfilename);
end

for i = 1:length(path_string)

    if ~exist(path_string{i}, 'dir')
        error('bds:PathNotExist', 'The string %s does not correspond to an existing directory.', path_string);
    end
    
    addpath(path_string{i});
    
    % Try saving the path in the system path-defining file at sys_pathdef. If the user does not have
    % writing permission for this file, then the path will not saved.
    % N.B. Do not save the path to the pathdef.m file under userpath. This file is not loaded by default
    % at startup. See
    % https://www.mathworks.com/matlabcentral/answers/269482-is-userdata-pathdef-m-for-local-path-additions-supported-on-linux
    orig_warning_state = warning;
    warning('off', 'MATLAB:SavePath:PathNotSaved'); % Maybe we do not have the permission to save path
    sys_pathdef = fullfile(matlabroot(), 'toolbox', 'local', 'pathdef.m');
    path_saved = (savepath(sys_pathdef) == 0);
    warning(orig_warning_state); % Restore the behavior of displaying warnings
    
    % If path not saved, try editing the startup.m of this user. Do this only if userpath is nonempty.
    % Otherwise, we will only get a startup.m in the current directory, which will not be executed
    % when MATLAB starts from other directories. On Linux, the default value of userpath is
    % ~/Documents/MATLAB, but it will be '' if this directory does not exist. We refrain from creating
    % this directory in that case.
    if ~path_saved && numel(userpath) > 0
        user_startup = fullfile(userpath, 'startup.m');
        add_path_string = sprintf('addpath(''%s'');', path_string{i});
        full_add_path_string = sprintf('%s  %s %s', add_path_string, '%', path_string_stamp);
    
        % First, check whether full_add_path_string already exists in user_startup or not.
        if exist(user_startup, 'file')
            startup_text_cells = regexp(fileread(user_startup), '\n', 'split');
            path_saved = any(strcmp(startup_text_cells, full_add_path_string));
        end
    
        if ~path_saved
            % We first check whether the last line of the user startup script is an empty line (or the
            % file is empty or even does not exist at all). If yes, we do not need to put a line break
            % before the path adding string.
            if exist(user_startup, 'file')
                startup_text_cells = regexp(fileread(user_startup), '\n', 'split');
                last_line_empty = isempty(startup_text_cells) || (isempty(startup_text_cells{end}) && ...
                    isempty(startup_text_cells{max(1, end-1)}));
            else
                last_line_empty = true;
            end
    
            file_id = fopen(user_startup, 'a');  % Open/create file for writing. Append data to the end.
            if file_id ~= -1 % If FOPEN cannot open the file, it returns -1; We keep silent if it fails.
                if ~last_line_empty  % The last line of user_startup is not empty
                    fprintf(file_id, '\n');  % Add a new empty line
                end
                fprintf(file_id, '%s', full_add_path_string);
                fclose(file_id);
                % Check that full_add_path_string is indeed added to user_startup.
                if exist(user_startup, 'file')
                    startup_text_cells = regexp(fileread(user_startup), '\n', 'split');
                    path_saved = any(strcmp(startup_text_cells, full_add_path_string));
                end
            end
        end
    end
end


% add_save_path ends
return
