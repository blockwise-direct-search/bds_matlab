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

fullpath = mfilename("fullpath");
path_tests = fileparts(fullpath);
path_bds = fileparts(path_tests);
path_src = fullfile(path_bds, "src");
path_competitors = fullfile(path_bds, "tests", "competitors");
addpath(path_tests)
addpath(path_src)
addpath(path_competitors)

cd(path_tests)
cd ./private
locate_matcutest();
cd (path_tests)

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

test_options = struct();
if ~isfield(parameters, "noise_type")
    test_options.noise_type = "gaussian";
end

if isfield(parameters, "feature")
    if startsWith(lower(parameters.feature), "randomx0")
        test_options.is_noisy = false;
        test_options.random_initial_point = true;
        level_str = split(lower(parameters.feature), "_");
        parameters.x0_perturbation_level = str2double(level_str{2});
    else
        switch lower(parameters.feature)
            case "plain"
                test_options.is_noisy = false;
                test_options.noise_level = 0;
                parameters.feature = "no_noise";
            case "negligible"
                test_options.is_noisy = true;
                test_options.noise_level = 1.0e-7;
                parameters.feature = strcat(test_options.noise_type, "_", "-7", "_noise");
            case "low"
                test_options.is_noisy = true;
                test_options.noise_level = 1.0e-5;
                parameters.feature = strcat(test_options.noise_type, "_", "-5", "_noise");
            case "medium"
                test_options.is_noisy = true;
                test_options.noise_level = 1.0e-3;
                parameters.feature = strcat(test_options.noise_type, "_", "-3", "_noise");
            case "high"
                test_options.is_noisy = true;
                test_options.noise_level = 1.0e-1;
                parameters.feature = strcat(test_options.noise_type, "_", "-1", "_noise");
            otherwise
                error("Unknown feature %s", parameters.feature);
        end
    end
end

if ~isfield(test_options, 'is_noisy')
    test_options.is_noisy = false;
end

if ~isfield(test_options, 'is_abs_noise')
    test_options.is_abs_noise = false;
end

if ~isfield(parameters, "log_x_axis")
    parameters.log_x_axis = false;
end

if ~isfield(parameters, "feature")
    parameters.feature = '';
end

if parameters.log_x_axis
    stamp = strcat(pdfname, "_", "mindim", "_", num2str(mindim), "_", ...
        "maxdim", "_", num2str(maxdim), "_", parameters.feature, "_", "log_x_axis", "_", time_str);
else
    stamp = strcat(pdfname, "_", "mindim", "_", num2str(mindim), "_", ...
        "maxdim", "_", num2str(maxdim), "_", parameters.feature, "_", time_str);
end

savepath = fullfile(path_testdata, stamp);
mkdir(savepath);
parameters.savepath = savepath;

prob_list = dimensions(type, mindim, maxdim);

for i = 1:length(prob_list)
    problem_name = prob_list{i};
    fhist_solvers(problem_name, parameters, test_options);
end

outdir = savepath;
outputfile = char(strcat("merged", "_", stamp, ".pdf"));
compdf_location = char(fullfile(path_tests, "private", "compdf"));
merge_pdf(outdir, outputfile, compdf_location);

cd(path_tests);
rmpath(path_tests);
rmpath(path_src);
rmpath(path_competitors);

end

