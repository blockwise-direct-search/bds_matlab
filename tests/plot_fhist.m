function plot_fhist(parameters)


if nargin < 3
    parameters = struct();
end

if isfield(parameters, "type")
    type = parameters.type;
else
    type = 'u';
end

if isfield(parameters, "mindim")
    mindim = parameters.mindim;
else
    mindim = 10;
end

if isfield(parameters, "maxdim")
    maxdim = parameters.maxdim;
else
    maxdim = 10;
end

% If the folder of testdata does not exist, make a new one.
current_path = mfilename("fullpath");
path_tests = fileparts(current_path);
path_testdata = fullfile(path_tests, "testdata");
if ~exist(path_testdata, "dir")
    mkdir(path_testdata);
end

% Use time to distinguish.
time_str = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm'));
% Trim time string.
time_str = trim_time(time_str);
stamp = strcat("mindim", "_", num2str(mindim), "_", "maxdim", "_", num2str(maxdim), "_", time_str);
savepath = fullfile(path_testdata, stamp);
mkdir(savepath);

parameters.solvers_name = ["bds", "newuoa"];
parameters.savepath = savepath;

prob_list = dimensions(type, mindim, maxdim);

for i = 1:length(prob_list)
    problem_name = prob_list{i};
    fhist_solvers(problem_name, parameters);
end

end

