function [best_value] = multi_tune(initial_value, parameters, options)
% This file is to tune the hyperparameters of BDS, including expanding
% factor, shrinking factor, reduction factors, shuffling_period and the
% replacement_delay.
%
% Input:
% Parameters should include the following fields:
% solvers_name: the algorithm of the BDS, including "cbds", "pbds", "rbds".
% problem_mindim: the minimum dimension of the problem. The default value
%                 is 1.
% problem_maxdim: the maximum dimension of the problem. The default value
%                 is 5.
% tau: the tolerance for performance profile.
% min_precision: the minimum precision for performance profile.
% parallel: whether to use parallel computing. The default value is false.
%
% Options should include the following fields:
% maxfun: the maximum number of function evaluations.
% initial_value: the initial value of the hyperparameters. If the algorithm
%                is cbds, the initial value should be a vector of length 5.
%                If the algorithm is pbds or rbds, the initial value should
%                be a vector of length 6 (including the shuffling_period
%                or replacement_delay).
%

% Record the current path.
oldpath = path();
% Restore the "right out of the box" path of MATLAB.
restoredefaultpath;
% Record the current directory.
old_dir = pwd();

% Add the paths that we need to use in the performance profile into the MATLAB
% search path.
current_path = mfilename("fullpath");
path_tests = fileparts(current_path);
path_root = fileparts(path_tests);
path_src = fullfile(path_root, "src");
path_competitors = fullfile(path_tests, "competitors");
addpath(path_root);
addpath(path_tests);
addpath(path_src);
addpath(path_competitors);

% If the folder of testdata does not exist, make a new one.
path_testdata = fullfile(path_tests, "testdata");
if ~exist(path_testdata, "dir")
    mkdir(path_testdata);
end

% We put time_str and tst here for the sake of plot_fhist.
% Use time to distinguish.
time_str = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm'));
% Trim time string.
time_str = trim_time(time_str);
tst = strcat(time_str, "_", "hyperparameters_tuning");

% Setup prima.
locate_prima_position = fullfile(path_competitors, "private");
cd(locate_prima_position)
locate_prima();

% Setup the constraints.
% The first index in lb and ub corresponds to the expanding factor.
% The second index in lb and ub corresponds to the shrinking factor.
% The values corresponding to indices 3 to 5 point to reduction_factors.
% We already know that reduction factors should not be large, so we set
% the upper bound as 1.
lb = [1, eps, 0, eps, eps];
ub = [10, 1-eps, 1, 1, 1];
Aineq = [0, 0, 1, -1, 0; 0, 0, 0, 1, -1];
bineq = [0; 0];
Aeq = [];
beq = [];

% Go back to the original path.
cd(old_dir);

% Set tau for performance profile.
if isfield(parameters, "min_precision")
    parameters.tau = 10.^(-1:-1:(-parameters.min_precision));
elseif ~isfield(parameters, "min_precision") && ~isfield(parameters, "tau")
    parameters.tau = 10.^(-1:-1:get_default_profile_options("min_precision"));
end
if ~isfield(parameters, "parallel")
    parameters.parallel = false;
end

options.output_xhist = true;

if strcmpi(parameters.solvers_name(1), "cbds")
    best_value = NaN(1, 7);
    % Preconditions for lincoa.
    initial_value = log(initial_value);
end
best_value(1) = length(parameters.tau);

if strcmpi(parameters.solvers_name(1), "cbds")
    % Here we should use lincoa since the constraint for the hyperparameters
    % is linear, including the expanding factor, shrinking factor, and the reduction
    % factors.
    [xopt, fopt, ~, output_tuning] = ...
        lincoa(@(x)perfprof_handle(exp(x), parameters), initial_value, Aineq, bineq, Aeq, beq, log(lb), log(ub), options);
    best_value(2) = fopt;
    best_value(3:7) = exp(xopt');
    % Scale xhist as the real value.
    output_tuning.xhist = exp(output_tuning.xhist);
end

path_testdata = fullfile(path_tests, "testdata");
path_testdata_outdir = fullfile(path_tests, "testdata", tst);

% Make a new folder to save numerical results and source code.
mkdir(path_testdata, tst);
mkdir(path_testdata_outdir, "tuning_data");
path_testdata_tuning_data = fullfile(path_testdata_outdir, "tuning_data");
mkdir(path_testdata_outdir, "src");
path_testdata_src = fullfile(path_testdata_outdir, "src");
mkdir(path_testdata_outdir, "tests");
path_testdata_tests = fullfile(path_testdata_outdir, "tests");
path_testdata_competitors = fullfile(path_testdata_tests, "competitors");
mkdir(path_testdata_competitors);
path_testdata_private = fullfile(path_testdata_tests, "private");
mkdir(path_testdata_private);

% Make a Txt file to store the initial_value, best_value and output under
% different tau.
filePath = strcat(path_testdata_tuning_data, "/tune_results.txt");
fileID = fopen(filePath, 'w');
% Write field names and their corresponding values into a file line by line.

best_value_record = num2str(best_value(i_tau, :));
separator = ", ";
best_value_record = strjoin(strsplit(best_value_record), separator);
fprintf(fileID, '%s: %s\n',"tau_length", best_value_record);
% Get the field names of the output structure under the specific tau.
output_tuning_saved = output_tuning{i_tau};
% Normally, output_tuning_saved should be a structure, which contains
% the following fields:
% xhist: the history of the hyperparameters.
% fhist: the history of the objective function values, which are the
%        value of the performance profile.
% funcCount: the number of function evaluations.
% algorithm: the algorithm used.
% message: the message returned by the algorithm.
% After trimming, the fields of output_tuning_saved should be still
% the same. However, the values of the fields may be different.
% xhist will be unchanged since it is always the matrix (here it implies
% that the number of rows are not 1).
% fhist will become a string. We will write it to the txt file
% directly.
% funcCount will also become a string. We will write it to the txt file
% directly, too.
% The same work for algorithm and message.
output_tuning_saved = trim_struct(output_tuning_saved);
% Get the field names of a structure.
output_tuning_saved_fields = fieldnames(output_tuning_saved);
for i = 1:numel(output_tuning_saved_fields)
    field = output_tuning_saved_fields{i};
    value = output_tuning_saved.(field);
    if ~iscell(value)
        if isnumvec(value)
            if size(value{1}, 1) == 1
                value_saved = num2str(value);
                separator = ", ";
                if length(value_saved) ~= 1
                    value_saved = strjoin(strsplit(value_saved), separator);
                    fprintf(fileID, '%s: %s\n', field, value_saved);
                end
            end
            % This extreme case is just for xhist, since it is still a matrix
            % after trimming.
        elseif ismatrix(value) && size(value, 1) == length(initial_value)
            fprintf(fileID, '%s:\n', field);
            % Define the width of each column and the spacing between columns.
            column_width = 10;
            column_spacing = 5;
            [rows, cols] = size(value);
            for row = 1:rows
                for col = 1:cols
                    % Use a formatting string to set left alignment and
                    % spacing for columns.
                    format = sprintf('%%-%ds', column_width);
                    fprintf(fileID, format, num2str(value(row, col)));
                    % Add spacing between columns.
                    if col < cols
                        fprintf(fileID, repmat(' ', 1, column_spacing));
                    end
                end
                fprintf(fileID, '\n');
            end
        elseif ischarstr(value)
            fprintf(fileID, '%s: %s\n', field, value);
        end
    else
        fprintf(fileID, '%s:\n', field);
        for j = 1:length(value)
            fprintf(fileID, '%s\n', value{j});
        end
    end
end
fclose(fileID);

% Make a Txt file to store the parameters that are used.
filePath = strcat(path_testdata_tuning_data, "/parameters.txt");
fileID = fopen(filePath, 'w');
parameters_saved = parameters;
parameters_saved = trim_struct(parameters_saved);
% Get the field names of a structure.
parameters_saved_fields = fieldnames(parameters_saved);
% Write field names and their corresponding values into a file line by line.
for i = 1:numel(parameters_saved_fields)
    field = parameters_saved_fields{i};
    value = parameters_saved.(field);
    if ~iscell(value)
        fprintf(fileID, '%s: %s\n', field, value);
    else
        for j = 1:length(value)
            solvers_options_saved = trim_struct(value{j});
            solvers_options_saved_fields = fieldnames(solvers_options_saved);
            for k = 1:numel(solvers_options_saved_fields)
                solvers_options_saved_field = solvers_options_saved_fields{k};
                solvers_options_saved_value = solvers_options_saved.(solvers_options_saved_field);
                fprintf(fileID, '%s: %s ', solvers_options_saved_field, ...
                    solvers_options_saved_value);
            end
            fprintf(fileID, '\n');
        end
    end
end
fclose(fileID);

% Copy the source code and test code to path_outdir.
copyfile(fullfile(path_src, "*"), path_testdata_src);
copyfile(fullfile(path_competitors, "*"), path_testdata_competitors);
copyfile(fullfile(path_tests, "private", "*"), path_testdata_private);
copyfile(fullfile(path_root, "setup.m"), path_testdata_outdir);

source_folder = path_tests;
destination_folder = path_testdata_tests;

% Get all files in the source folder.
file_list = dir(fullfile(source_folder, '*.*'));
file_list = file_list(~[file_list.isdir]);

% Copy all files (excluding subfolders) to the destination folder.
for i = 1:numel(file_list)
    source_file = fullfile(source_folder, file_list(i).name);
    destination_file = fullfile(destination_folder, file_list(i).name);
    copyfile(source_file, destination_file);
end

% Plot the performance profile with the specific tau.
parameters_perfprof.solvers_name = parameters.solvers_name;
parameters_perfprof.parallel = true;
if ~isfield(parameters, "problem_mindim")
    parameters_perfprof.problem_mindim = 1;
else
    parameters_perfprof.problem_mindim = parameters.problem_mindim;
end
if ~isfield(parameters, "problem_maxdim")
    parameters_perfprof.problem_maxdim = 5;
else
    parameters_perfprof.problem_maxdim = parameters.problem_maxdim;
end
for i = 1:length(tau)
    parameters_perfprof.tau = tau_tuning(i);
    if strcmpi(parameters_perfprof.solvers_name(1), "cbds")
        parameters_perfprof.solvers_options{1}.expand = best_value(i, 3);
        parameters_perfprof.solvers_options{1}.shrink = best_value(i, 4);
        parameters_perfprof.solvers_options{1}.reduction_factor = best_value(i, 5:7);
        plot_profile(parameters_perfprof);
    end
end

% Restore the path to oldpath.
setpath(oldpath);
% Go back to the original directory.
cd(old_dir);

end

