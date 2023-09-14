function [output] = profile(parameters)
% Draw performance profiles.

% Record the current path.
oldpath = path();
% Restore the "right out of the box" path of MATLAB.
restoredefaultpath;
% Record the current directory.
old_dir = pwd();

exception = [];

try

    path_nlopt = "/usr/local/lib/matlab/";

    if ~contains(path, path_nlopt, 'IgnoreCase', true)
        if exist(path_nlopt, 'dir') == 7
            addpath(path_nlopt);
            disp('Add path_nlopt to MATLAB path.');
        else
            disp('path_nlopt does not exist on the local machine.');
        end

    else
        disp('Path_nlopt already exists in MATLAB path.');
    end

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

    % In case no solvers are input, then throw an error.
    if ~isfield(parameters, "solvers_options") || length(parameters.solvers_options) < 2
        error("There should be at least two solvers.")
    end

    % Get the parameters that the test needs.
    parameters = set_profile_options(parameters);

    % Tell MATLAB where to find MatCUTEst.
    locate_matcutest();
    % Tell MATLAB where to find prima.
    locate_prima();

    % Get list of problems
    s.type = parameters.problems_type; % Unconstrained: 'u'
    s.mindim = parameters.problems_mindim; % Minimum of dimension
    s.maxdim = parameters.problems_maxdim; % Maximum of dimension
    s.blacklist = [];
    %s.blacklist = [s.blacklist, {}];
    % Problems that takes too long to solve.
    solvers_name = [];
    for i = 1:length(parameters.solvers_options)
        solvers_name = [solvers_name, parameters.solvers_options{i}.solver];
    end
    
    % {'FBRAIN3LS'} and {'STRATEC'} take too long for fminunc.
    if ismember("matlab_fminunc", solvers_name)
        s.blacklist = [s.blacklist, {'FBRAIN3LS'}, {'STRATEC'}];
    end
    % {"MUONSINELS"} takes nlopt_newuoa so long to run (even making MATLAB crash).
    % {"LRCOVTYPE"}, {'HIMMELBH'} and {'HAIRY'} take nlopt_cobyla so long
    % to run (even making MATLAB crash).
    % {"MUONSINELS"} takes nlopt_bobyqa so long to run (even making MATLAB crash).
    if ismember("nlopt", solvers_name)
        s.blacklist = [s.blacklist, {'MUONSINELS'}, {'BENNETT5LS'},...
            {'HIMMELBH'}, {'HAIRY'}];
    end
    
    problem_names = secup(s);

    fprintf("We will load %d problems\n\n", length(problem_names))

    % Some fixed (relatively) options
    % TODO: Read two papers: What Every Computer Scientist Should Know About
    % Floating-Point Arithmetic; stability and accuracy numerical(written by Higham).

    % Set maxfun for frec.
    if isfield(parameters, "maxfun_factor") && isfield(parameters, "maxfun")
        maxfun_frec = max(parameters.maxfun_factor*parameters.problems_maxdim, parameters.maxfun);
    elseif isfield(parameters, "maxfun_factor")
        maxfun_frec = parameters.maxfun_factor*parameters.problems_maxdim;
    elseif isfield(parameters, "maxfun")
        maxfun_frec = parameters.maxfun;
    else
        maxfun_frec = max(get_default_profile_options("maxfun"), ...
            get_default_profile_options("maxfun_factor")*parameters.problems_maxdim);
    end

    % Initialize fmin and frec.
    num_solvers = length(parameters.solvers_options);
    % Get number of problems.
    num_problems = length(problem_names);
    % Get Number of random tests(If num_random = 1, it means no random test).
    num_random = parameters.num_random;
    % Record minimum value of the problems of the random test.
    fmin = NaN(num_problems, num_random);
    frec = NaN(num_problems, num_solvers, num_random, maxfun_frec);

    % Set noisy parts of test.
    test_options.is_noisy = parameters.is_noisy;
    if parameters.is_noisy

        if isfield(parameters, "noise_level")
            test_options.noise_level = parameters.noise_level;
        else
            test_options.noise_level = get_default_profile_options("noise_level");
        end

        % Relative: (1+noise_level*noise)*f; absolute: f+noise_level*noise
        if isfield(parameters, "is_abs_noise")
            test_options.is_abs_noise = parameters.is_abs_noise;
        else
            test_options.is_abs_noise = get_default_profile_options("is_abs_noise");
        end

        if isfield(parameters, "noise_type")
            test_options.noise_type = parameters.noise_type;
        else
            test_options.noise_type = get_default_profile_options("noise_type");
        end

        if isfield(parameters, "num_random")
            test_options.num_random = parameters.num_random;
        else
            test_options.num_random = get_default_profile_options("num_random");
        end

    end
    
    % Set maxfun and maxfun_factor for test_options.
    if isfield(parameters, "maxfun")
        test_options.maxfun = parameters.maxfun;
    end

    if isfield(parameters, "maxfun_factor")
        test_options.maxfun_factor = parameters.maxfun_factor;
    end

    % Set scaling matrix.
    test_options.scale_variable = false;

    % Set solvers_options.
    solvers_options = parameters.solvers_options;

    % If parameters.noise_initial_point is true, then initial point will be
    % selected for each problem num_random times.
    % parameters.fmintype is set to be "randomized" defaultly, then there is
    % no need to test without noise, which makes the curve of performance profile
    % more higher. If parallel is true, use parfor to calculate (parallel computation),
    % otherwise, use for to calculate (sequential computation).
    if parameters.parallel == true
        parfor i_problem = 1:num_problems
            p = macup(problem_names(1, i_problem));
            for i_run = 1:num_random
                fval_tmp = NaN(1, num_solvers);
                if parameters.random_initial_point
                    rr = randn(size(p.x0));
                    rr = rr / norm(rr);
                    p.x0 = p.x0 + 1e-3 * max(1, norm(p.x0)) * rr;
                end
                fprintf("%d(%d). %s\n", i_problem, i_run, p.name);
                for i_solver = 1:num_solvers
                    fhist = get_fhist(p, maxfun_frec, i_solver,...
                        i_run, solvers_options, test_options);
                    fval_tmp(i_solver) = min(fhist);
                    frec(i_problem,i_solver,i_run,:) = fhist;
                end
                fmin(i_problem, i_run) = min(fval_tmp);
            end
        end
    else
        for i_problem = 1:num_problems
            p = macup(problem_names(1, i_problem));
            for i_run = 1:num_random
                fval_tmp = NaN(1, num_solvers);
                if parameters.random_initial_point
                    rr = randn(size(p.x0));
                    rr = rr / norm(rr);
                    p.x0 = p.x0 + 1e-3 * max(1, norm(p.x0)) * rr;
                end
                fprintf("%d(%d). %s\n", i_problem, i_run, p.name);
                for i_solver = 1:num_solvers
                    fhist = get_fhist(p, maxfun_frec, i_solver,...
                        i_run, solvers_options, test_options);
                    fval_tmp(i_solver) = min(fhist);
                    frec(i_problem,i_solver,i_run,:) = fhist;
                end
                fmin(i_problem, i_run) = min(fval_tmp);
            end
        end
    end

    % If parameters.fmintype = "real-randomized", then test without noise
    % should be conducted and fmin might be smaller, which makes curves
    %  of performance profile more lower.
    if test_options.is_noisy && strcmpi(parameters.fmin_type, "real-randomized")
        fmin_real = NaN(num_problems, 1);
        test_options.is_noisy = false;
        i_run = 1;
        if parameters.parallel == true
            parfor i_problem = 1:num_problems
                p = macup(problem_names(1, i_problem));
                frec_local = NaN(num_solvers, maxfun_frec);
                if parameters.random_initial_point
                    rr = randn(size(x0));
                    rr = rr / norm(rr);
                    p.x0 = p.x0 + 1e-5 * max(1, norm(p.x0)) * rr;
                end
                fprintf("%d. %s\n", i_problem, p.name);
                for i_solver = 1:num_solvers
                    frec_local(i_solver,:) = get_fhist(p, maxfun_frec,...
                        i_solver, i_run, solvers_options, test_options);
                end
                fmin_real(i_problem) = min(frec_local(:, :),[],"all");
            end
        else
            for i_problem = 1:num_problems
                p = macup(problem_names(1, i_problem));
                frec_local = NaN(num_solvers, maxfun_frec);
                if parameters.random_initial_point
                    rr = randn(size(x0));
                    rr = rr / norm(rr);
                    p.x0 = p.x0 + 1e-3 * max(1, norm(p.x0)) * rr;
                end
                fprintf("%d. %s\n", i_problem, p.name);
                for i_solver = 1:num_solvers
                    frec_local(i_solver,:) = get_fhist(p, maxfun_frec,...
                        i_solver, i_run, solvers_options, test_options);
                end
                fmin_real(i_problem) = min(frec_local(:, :),[],"all");
            end
        end
    end

    if strcmpi(parameters.fmin_type, "real-randomized")
        fmin_total = [fmin, fmin_real];
        fmin = min(fmin_total, [], 2);
    end

    % Use time to distinguish.
    time = datetime("now");
    time_str = sprintf('%04d-%02d-%02d %02d:%02d:%02.0f', year(time), ...
        month(time), day(time), hour(time), minute(time), second(time));
    % Trim time string.
    time_str = trim_time(time_str);
    tst = sprintf("test_%s", time_str);
    % Rename tst as mixture of time stamp and pdfname.
    tst = strcat(tst, "_", parameters.pdfname);
    path_testdata = fullfile(path_tests, "testdata");
    path_testdata_outdir = fullfile(path_tests, "testdata", tst);

    % Make a new folder to save numerical results and source code.
    mkdir(path_testdata, tst);
    mkdir(path_testdata_outdir, "perf");
    options_perf.outdir = fullfile(path_testdata_outdir, "perf");
    mkdir(path_testdata_outdir, "src");
    path_testdata_src = fullfile(path_testdata_outdir, "src");
    mkdir(path_testdata_outdir, "tests");
    path_testdata_tests = fullfile(path_testdata_outdir, "tests");
    path_testdata_competitors = fullfile(path_testdata_tests, "competitors");
    mkdir(path_testdata_competitors);
    path_testdata_private = fullfile(path_testdata_tests, "private");
    mkdir(path_testdata_private);

    % Copy the source code and test code to path_outdir.
    copyfile(fullfile(path_src, "*"), path_testdata_src);
    copyfile(fullfile(path_competitors, "*"), path_testdata_competitors);
    copyfile(fullfile(path_tests, "private", "*"), path_testdata_private);

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

    % Draw performance profiles.
    % Set tolerance of convergence test in performance profile.
    tau = parameters.tau;
    tau_length = length(tau);
    options_perf.time_stamp = time_str;
    options_perf.solvers = parameters.solvers_legend;
    options_perf.natural_stop = false;

    for l = 1:tau_length
        options_perf.tau = tau(l);
        output = perfprof(frec, fmin, options_perf);
    end

    cd(options_perf.outdir);

    % Initialize string variable.
    pdfFiles = dir(fullfile(options_perf.outdir, '*.pdf'));

    % Store filename in a cell.
    pdfNamesCell = cell(numel(pdfFiles), 1);
    for i = 1:numel(pdfFiles)
        pdfNamesCell{i} = pdfFiles(i).name;
    end

    % Use the strjoin function to concatenate the elements in a cell array into a single string.
    inputfiles = strjoin(pdfNamesCell, ' ');

    % Remove spaces at the beginning of a string.
    inputfiles = strtrim(inputfiles);

    % Merge pdf.
    outputfile = 'all.pdf';
    system(['bash ', fullfile(path_tests, 'private', 'compdf'), ' ', inputfiles, ' -o ', outputfile]);
    % Rename pdf.
    movefile("all.pdf", sprintf("%s.pdf", parameters.pdfname));

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
