function plot_fhist(parameters)

% In case no solvers are input, then throw an error.
if ~isfield(parameters, "solvers_name") || length(parameters.solvers_name) < 2
    error("There should be at least two solvers.")
end

for i = 1:length(parameters.solvers_name)
    parameters.solvers_options{i}.solver = parameters.solvers_name(i);
end

parameters = get_solvers(parameters);
%parameters = rmfield(parameters, "solvers_name");

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
for i = 1:length(parameters.solvers_name)
    pdfname_solver = get_pdf_name(parameters, i);
    if i == 1
        pdfname = pdfname_solver;
    else
        pdfname = strcat(pdfname, "_", pdfname_solver);
    end
end

stamp = strcat(pdfname, "_", "mindim", "_", num2str(mindim), "_", "maxdim", "_", num2str(maxdim), "_", time_str);
savepath = fullfile(path_testdata, stamp);
mkdir(savepath);
parameters.savepath = savepath;

prob_list = dimensions(type, mindim, maxdim);

for i = 1:length(prob_list)
    problem_name = prob_list{i};
    fhist_solvers(problem_name, parameters);
end

end

