function verify(parameters)

% Record the current path.
oldpath = path();
% Restore the "right out of the box" path of MATLAB.
restoredefaultpath;
% Record the current directory.
old_dir = pwd();

exception = [];

try

    % In case no version of archiva is input.
    if ~isfield(parameters, "version")
        error("The version of archiva should be input.")
    end

    % Compile the verison of archiva.
    path_norma = locate_norma();
    path_version = fullfile(path_norma, parameters.version);
    path_verify = fileparts(mfilename('fullpath'));
    cd(path_version);
    setup
    cd(path_verify);

    % Compile the version of modern repository.
    path_root = fileparts(path_verify);
    cd(path_root);
    setup
    cd(path_verify);

    solvers = {"bds", "bds_norma"};

    % Tell MATLAB where to find MatCUTEst.
    locate_matcutest();

    % Get list of problems
    if isfield(parameters, "problems_type")
        s.type = parameters.problems_type;
    else
        s.type = get_default_profile_options("problems_type");
    end

    if isfield(parameters, "problems_mindim")
        s.mindim = parameters.problems_mindim;
    else
        s.mindim = get_default_profile_options("problems_mindim");
    end

    if isfield(parameters, "problems_maxdim")
        s.maxdim = parameters.problems_maxdim;
    else
        s.maxdim = 20;
    end
    
    s.blacklist = [];

    if s.mindim >= 6
        s.blacklist = [parameters.blacklist, { 'ARGTRIGLS', 'BROWNAL', ...
            'COATING', 'DIAMON2DLS', 'DIAMON3DLS', 'DMN15102LS', ...
            'DMN15103LS', 'DMN15332LS', 'DMN15333LS', 'DMN37142LS', ...
            'DMN37143LS', 'ERRINRSM', 'HYDC20LS', 'LRA9A', ...
            'LRCOVTYPE', 'LUKSAN12LS', 'LUKSAN14LS', 'LUKSAN17LS', 'LUKSAN21LS', ...
            'LUKSAN22LS', 'MANCINO', 'PENALTY2', 'PENALTY3', 'VARDIM',
            }];
    end

    problem_names = secup(s);

    fprintf("We will load %d problems\n\n", length(problem_names))

    % Get the number of problems.
    num_problems = length(problem_names);

    if isfield(parameters, 'prec')
        prec = parameters.prec;
    else
        prec = 0;
    end

    if isfield(parameters, 'single_test')
        single_test = parameters.single_test;
    else
        single_test = false;
    end

    if ~isfield(parameters, "parallel")
        parallel = get_default_profile_options("parallel");
    end

    if ~isfield(parameters, "num_random")
        num_random = get_default_profile_options("num_random");
    end

    if parallel
        parfor i_problem = 1:num_problems
            p = macup(problem_names(1, i_problem));
            for i_run = 1:num_random
                fprintf("%d(%d). %s\n", i_problem, i_run, p.name);
                iseqiv(solvers, p, i_run, prec, single_test, parameters);
            end
        end
    else
        for i_problem = 1:num_problems
            p = macup(problem_names(1, i_problem));
            for i_run = 1:num_random
                fprintf("%d(%d). %s\n", i_problem, i_run, p.name);
                iseqiv(solvers, p, i_run, prec, single_test, parameters);
            end
        end
    end

    % Use time to distinguish.
    time_str = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm'));
    % Trim time string.
    time_str = trim_time(time_str);
    % tst = sprintf("test_%s", time_str);
    % Rename tst as the mixture of time stamp and pdfname.
    tst = strcat(parameters.pdfname, "_", time_str);
    path_testdata = fullfile(path_tests, "testdata");
    path_testdata_outdir = fullfile(path_tests, "testdata", tst);

    % Make a new folder to save numerical results and source code.
    mkdir(path_testdata, tst);
    mkdir(path_testdata_outdir, "src");
    path_testdata_src = fullfile(path_testdata_outdir, "src");

    % Make a Txt file to store the parameters that are used.
    filePath = strcat(path_testdata_perf, "/parameters.txt");
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

catch exception

    % Do nothing for the moment.

end

% Restore the path to oldpath.
setpath(oldpath);
% Go back to the original directory.
cd(old_dir);

if ~isempty(exception)  % Rethrow any exception caught above.
    rethrow(exception);
end

end