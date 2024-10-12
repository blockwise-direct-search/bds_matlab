function profile_multi_feature(parameters)
% In profile_multi_feature, the required field of parameters are:
%  solver_name: an array of strings, the names of the solvers to be profiled, like ["cbds", "newuoa"].
%  multi_feature: an array of strings, the names of the features to be profiled, like ["plain", "rotation"].
%  If it is "all", all features will be profiled and the features going to be profiled are the following:
%  ["plain", "rotation", "randomx0_10", "rotation_noisy_1e-1", "rotation_noisy_1e-2", "rotation_noisy_1e-3", "rotation_noisy_1e-4", "rotation_noisy_1e-5"].
%  problem_dim: a string, the dimension of the problems to be profiled, like "small" or "big".
%

multi_feature = parameters.multi_feature;
if strcmpi(parameters.multi_feature, "all")
    multi_feature = ["plain", "rotation", "randomx0_10", "rotation_noisy_1e-1",...
        "rotation_noisy_1e-2", "rotation_noisy_1e-3", "rotation_noisy_1e-4",...
        "rotation_noisy_1e-5", "rotation_noisy_1e-6"];
end
if isfield(parameters, "solvers_name") && any(contains(parameters.solvers_name, "bfgs"))
    multi_feature = multi_feature(~contains(multi_feature, "randomx0"));
end

parameters = rmfield(parameters, "multi_feature");
pdfname_feature = "";
for i = 1:length(multi_feature)
    if i == 1
        pdfname_feature = multi_feature{i};
    else
        pdfname_feature_base = strrep(multi_feature{i}, '''', '"');
        pdfname_feature = strcat(pdfname_feature, "_", pdfname_feature_base);
    end
end

if isfield(parameters, "problem_dim")
    if strcmpi(parameters.problem_dim, "small")
        parameters.problem_mindim = 1;
        parameters.problem_maxdim = 5;
    elseif strcmpi(parameters.problem_dim, "big")
        parameters.problem_mindim = 6;
        parameters.problem_maxdim = 200;
    end
end

if ~isfield(parameters, "fmin_type")
    parameters.fmin_type = get_default_profile_options("fmin_type");
end

if ~isfield(parameters, "test_type")
    parameters.test_type = "matcutest";
end

pdfname_solver = "";
for i = 1:length(parameters.solvers_name)
    if i == 1
        pdfname_solver = parameters.solvers_name{1};
    else
        pdfname_solver_base = strrep(parameters.solvers_name{i}, '''', '"');
        pdfname_solver = strcat(pdfname_solver, "_", pdfname_solver_base);
    end
end

pdfname = strcat("merged", "_", pdfname_solver, "_", num2str(parameters.problem_mindim), "_",...
    num2str(parameters.problem_maxdim), "_", parameters.fmin_type, "_", pdfname_feature, "_", parameters.test_type);

% We put time_str and tst here for the sake of plot_fhist.
% Use time to distinguish.
time_str = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm'));
% Trim time string.
time_str = trim_time(time_str);
% Rename tst as the mixture of time stamp and pdfname.
tst = strcat(pdfname, "_", time_str);

current_path = mfilename("fullpath");
path_tests = fileparts(current_path);
path_testdata = fullfile(path_tests, "testdata");
multi_feature_outdir = fullfile(path_testdata, tst);
mkdir(multi_feature_outdir);

for i = 1:length(multi_feature)
    parameters.feature = strrep(multi_feature{i}, '''', '"');

    if ~isfield(parameters, "solvers_options")
        parameters.solvers_options = {};
    end

    for j = 1:length(parameters.solvers_name)
        parameters.solvers_options{j}.solver = parameters.solvers_name(j);
    end

    [path_testdata_perf, ~, ~] = profile(parameters);

    % Search for PDF files starting with 'merge' in the source folder.
    pdf_files = dir(fullfile(path_testdata_perf, 'merge*.pdf'));
    
    % Iterate through the found PDF files and copy them to the target folder.
    for k = 1:numel(pdf_files)
        % The full path to the build source and target files.
        source_file = fullfile(path_testdata_perf, pdf_files(k).name);
        target_file = fullfile(multi_feature_outdir, pdf_files(k).name);

        % Copy Files.
        copyfile(source_file, target_file);

        % Displays a message that the copy is successful.
        fprintf('Copied %s to %s\n', pdf_files(k).name, multi_feature_outdir);
    end

end

compdf_location = char(fullfile(path_tests, "private", "compdf"));
pdfname = char(pdfname);
merge_pdf_order(multi_feature_outdir, pdfname, compdf_location);

end

