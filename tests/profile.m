function [path_testdata_perf, frec, fmin] = profile(parameters)
% Draw performance profiles.
%

% Record the current path.
oldpath = path();
% Restore the "right out of the box" path of MATLAB.
restoredefaultpath;
% Clear the MATLAB search path.
rehash toolboxcache;
% Record the current directory.
old_dir = pwd();

exception = [];
path_testdata_perf = " ";

try

    if any(ismember(parameters.solvers_name, "nomad"))
        addpath '/home/lhtian97/Documents/nomad/build/release/lib'
    end

    % Add the paths that we need to use in the performance profile into the MATLAB
    % search path.
    current_path = mfilename("fullpath");
    path_tests = fileparts(current_path);
    path_root = fileparts(path_tests);
    path_src = fullfile(path_root, "src");
    path_competitors = fullfile(path_tests, "competitors");
    path_tools = fullfile(path_tests, "tools");
    addpath(path_root);
    addpath(path_tests);
    addpath(path_src);
    addpath(path_competitors);
    addpath(path_tools);

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

    % Get list of problems
    if isfield(parameters, "test_type") && strcmpi(parameters.test_type, "matcutest")

        s.type = parameters.problem_type; % Unconstrained: 'u'
        s.mindim = parameters.problem_mindim; % Minimum of dimension
        s.maxdim = parameters.problem_maxdim; % Maximum of dimension
        s.blacklist = [];
        %s.blacklist = [s.blacklist, {}];
        % Problems that take too long to solve.
        % {'FBRAIN3LS'} and {'STRATEC'} take too long for fminunc.
        if ismember("matlab_fminunc", parameters.solvers_name)
            s.blacklist = [s.blacklist, {'FBRAIN3LS'}, {'STRATEC'}];
        end
        % {"MUONSINELS"} takes nlopt_newuoa so long to run (even making MATLAB crash).
        % {"LRCOVTYPE"}, {'HIMMELBH'} and {'HAIRY'} take nlopt_cobyla so long
        % to run (even making MATLAB crash).
        % {"MUONSINELS"} takes nlopt_bobyqa so long to run (even making MATLAB crash).
        if ismember("nlopt", parameters.solvers_name)
            s.blacklist = [s.blacklist, {'MUONSINELS'}, {'BENNETT5LS'},...
                {'HIMMELBH'}, {'HAIRY'}];
        end

        %if s.mindim >= 6
        s.blacklist = [s.blacklist, { 'ARGTRIGLS', 'BROWNAL', ...
            'COATING', 'DIAMON2DLS', 'DIAMON3DLS', 'DMN15102LS', ...
            'DMN15103LS', 'DMN15332LS', 'DMN15333LS', 'DMN37142LS', ...
            'DMN37143LS', 'ERRINRSM', 'HYDC20LS', 'LRA9A', ...
            'LRCOVTYPE', 'LUKSAN12LS', 'LUKSAN14LS', 'LUKSAN17LS', 'LUKSAN21LS', ...
            'LUKSAN22LS', 'MANCINO', 'PENALTY2', 'PENALTY3', 'VARDIM',
            }];
        %end

        if isfield(parameters, "problem_names")
            problem_names = parameters.problem_names;
        else
            problem_names = secup(s);
        end
    else
        path_optiprofiler = strcat(fileparts(fileparts(fileparts(fileparts(mfilename('fullpath'))))), ...
            "/local", "/optiprofiler", "/matlab", "/optiprofiler");
        if ~exist(path_optiprofiler, "dir")
            error("The path of optiprofiler does not exist.");
        end
        path_optiprofiler_probinfo_s2mpj = strcat(path_optiprofiler, "/probinfo_s2mpj");
        path_optiprofiler_matlab_problems = strcat(path_optiprofiler_probinfo_s2mpj, "/matlab_problems");
        path_optiprofiler_src = strcat(path_optiprofiler, "/src");
        path_optiprofiler_tests = strcat(path_optiprofiler, "/tests");
        addpath(path_optiprofiler);
        addpath(path_optiprofiler_probinfo_s2mpj);
        addpath(path_optiprofiler_matlab_problems);
        addpath(path_optiprofiler_src);
        addpath(path_optiprofiler_tests);
        blacklist = ["DIAMON2DLS",...
            "DIAMON2D",...
            "DIAMON3DLS",...
            "DIAMON3D",...
            "DMN15102LS",...
            "DMN15102",...
            "DMN15103LS",...
            "DMN15103",...
            "DMN15332LS",...
            "DMN15332",...
            "DMN15333LS",...
            "DMN15333",...
            "DMN37142LS",...
            "DMN37142",...
            "DMN37143LS",...
            "DMN37143",...
            "ROSSIMP3_mp"];
        blacklist_time_consuming = ["BAmL1SPLS",...
            "FBRAIN3LS",...
            "GAUSS1LS",...
            "GAUSS2LS",...
            "GAUSS3LS",...
            "HYDC20LS",...
            "HYDCAR6LS",...
            "LUKSAN11LS",...
            "LUKSAN12LS",...
            "LUKSAN13LS",...
            "LUKSAN14LS",...
            "LUKSAN17LS",...
            "LUKSAN21LS",...
            "LUKSAN22LS",...
            "METHANB8LS",...
            "METHANL8LS",...
            "SPINLS",...
            "VESUVIALS",...
            "VESUVIOLS",...
            "VESUVIOULS",...
            "YATP1CLS"];
        options_s2mpj.problem_type = 'u';
        options_s2mpj.mindim = parameters.problem_mindim;
        options_s2mpj.maxdim = parameters.problem_maxdim;
        problem_names_orig = s_select(options_s2mpj);
        problem_names = [];
        for i = 1:length(problem_names_orig)
            if ~ismember(problem_names_orig(i), blacklist) && ...
                    ~ismember(problem_names_orig(i), blacklist_time_consuming)
                problem_names = [problem_names, problem_names_orig(i)];
            end
        end
    end

    fprintf("We will load %d problems\n\n", length(problem_names));

    % Some fixed (relatively) options
    % Read two papers: What Every Computer Scientist Should Know About
    % Floating-Point Arithmetic; stability and accuracy numerical(written by Higham).

    % Initialize the number of solvers.
    num_solvers = length(parameters.solvers_options);
    % Set MaxFunctionEvaluations_frec for performance profile.
    MaxFunctionEvaluations_frec = max(get_default_profile_options("MaxFunctionEvaluations"), ...
        get_default_profile_options("MaxFunctionEvaluations_dim_factor")*parameters.problem_maxdim);
    for i = 1:num_solvers
        if isfield(parameters, "default") && parameters.default
            switch parameters.solvers_options{i}.solver
                case "bfo_wrapper"
                    MaxFunctionEvaluations_frec = max(MaxFunctionEvaluations_frec, 5000*parameters.problem_maxdim);
                case "fminsearch_wrapper"
                    MaxFunctionEvaluations_frec = max(MaxFunctionEvaluations_frec, 200*parameters.problem_maxdim);
                case "fminunc_wrapper"
                    MaxFunctionEvaluations_frec = max(MaxFunctionEvaluations_frec, 100*parameters.problem_maxdim);
                case "prima_wrapper"
                    MaxFunctionEvaluations_frec = max(MaxFunctionEvaluations_frec, 500*parameters.problem_maxdim);
            end
        end
    end

    % Get the number of problems.
    num_problems = length(problem_names);
    % Get Number of random tests(If num_random = 1, it means no random test).
    num_random = parameters.num_random;
    fprintf("num_random = %d\n", num_random);
    % Record the minimum value of the problems of the random test.
    fmin = NaN(num_problems, num_random);
    frec = NaN(num_problems, num_solvers, num_random, MaxFunctionEvaluations_frec);

    % Set noisy parts of the test.
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

    % Set solvers_options.
    parameters = get_options(parameters);
    solvers_options = parameters.solvers_options;

    % We put time_str and tst here for the sake of plot_fhist.
    % Use time to distinguish.
    time_str = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm'));
    % Trim time string.
    time_str = trim_time(time_str);
    % Rename tst as the mixture of time stamp and pdfname.
    tst = strcat(parameters.pdfname, "_", time_str);


    if isfield(parameters, "plot_fhist") && parameters.plot_fhist
        if ~isfield(parameters, "log_x_axis")
            parameters.log_x_axis = false;
        end

        if parameters.log_x_axis
            parameters.stamp_fhist = strcat(parameters.pdfname, "_", "log_x_axis", "_", time_str, "_fhist");
        else
            parameters.stamp_fhist = strcat(parameters.pdfname, "_", time_str, "_fhist");
        end

        savepath = fullfile(path_testdata, parameters.stamp_fhist);
        mkdir(savepath);
        if num_random == 1
            parameters.savepath = savepath;
        else
            parameters.savepath = cell(1, num_random);
            for i = 1:num_random
                parameters.savepath{i} = fullfile(savepath, sprintf("random_%d", i));
                mkdir(parameters.savepath{i});
            end
        end
    end

    % If parameters.noise_initial_point is true, then the initial point will be
    % selected for each problem num_random times.
    % The default value of parameters.fmin_type is set to be "randomized", then there is
    % no need to test without noise, which makes the curve of the performance profile
    % more higher. If parallel is true, use parfor to calculate (parallel computation),
    % otherwise, use for to calculate (sequential computation).
    if parameters.parallel == true
        parfor i_problem = 1:num_problems
            if isfield(parameters, "test_type") && strcmpi(parameters.test_type, "matcutest")
                p = macup(problem_names(1, i_problem));
            else
                problem_info = s_load(char(problem_names(i_problem)));
                p = s2mpj_wrapper(problem_info);
            end
            dim = length(p.x0);
            for i_run = 1:num_random
                % Set scaling matrix.
                if isfield(parameters, "feature") && strcmpi(parameters.feature, "badly_scaled")
                    % Badly_scaled is a flag to indicate whether the problem is badly scaled.
                    seed = min(abs(ceil(1e5*sin(1e9*sum(num_problems)))) + ...
                        abs(ceil(1e4 * sin(1e7*num_random))) + 5000 * num_random, 2^32 - 1 -num_random) + i_run;
                    rng(seed)
                    if isfield(parameters, "badly_scaled_sigma")
                        scale_matrix = diag(2.^(parameters.badly_scaled_sigma*randn(dim, 1)));
                    else
                        scale_matrix = diag(2.^(dim*randn(dim, 1)));
                    end
                    h = @(x) scale_matrix * x;
                    p.objective = @(x) p.objective(h(x));
                    p.x0 = p.x0 ./ diag(scale_matrix);
                end
                if isfield(parameters, "feature") && strcmpi(parameters.feature, "rotation_badly_scaled")
                    % Rotation_badly_scaled is a flag to indicate whether the problem is rotated and badly scaled.
                    seed = min(abs(ceil(1e5*sin(1e9*sum(num_problems)))) + ...
                        abs(ceil(1e4 * sin(1e7*num_random))) + 5000 * num_random, 2^32 - 1 -num_random) + i_run;
                    rng(seed)
                    [Q,R] = qr(randn(dim));
                    Q(:, diag(R) < 0) = -Q(:, diag(R) < 0);
                    rotation_matrix = Q;
                    h = @(x) rotation_matrix * x;
                    p.objective = @(x) p.objective(h(x));
                    p.x0 = (qr(rotation_matrix) \ eye(dim)) * p.x0;
                    scale_matrix = diag(2.^(dim*randn(dim, 1)));
                    h = @(x) scale_matrix * x;
                    p.objective = @(x) p.objective(h(x));
                    p.x0 = p.x0 ./ diag(scale_matrix);
                end
                if isfield(parameters, "feature") && (strcmpi(parameters.feature, "rotation") || ...
                        strcmpi(parameters.feature, "rotation_noisy"))
                    % Rotation is a flag to indicate whether the problem is rotated.
                    seed = min(abs(ceil(1e5*sin(1e9*sum(num_problems)))) + ...
                        abs(ceil(1e4 * sin(1e7*num_random))) + 5000 * num_random, 2^32 - 1 -num_random) + i_run;
                    rng(seed)
                    [Q,R] = qr(randn(dim));
                    Q(:, diag(R) < 0) = -Q(:, diag(R) < 0);
                    rotation_matrix = Q;
                    h = @(x) rotation_matrix * x;
                    p.objective = @(x) p.objective(h(x));
                    p.x0 = (qr(rotation_matrix) \ eye(dim)) * p.x0;
                end
                if isfield(parameters, "feature") && strcmpi(parameters.feature, "structured")
                    % Structured is a flag to indicate whether the problem is added with l-p regularization term.
                    if isfield(parameters, "structured_factor")
                        h = @(x) parameters.structured_factor *  sum(abs(x).^ 1);
                    else
                        h = @(x) sum(abs(x).^ 1);
                    end
                    %h = @(x) parameters.structured_factor *  sum(abs(x).^ parameters.structured_norm)^(1/parameters.structured_norm);
                    p.objective = @(x) p.objective(x) + h(x);
                end
                if isfield(parameters, "feature") && strcmpi(parameters.feature, "rotation_structured")
                    % Rotation_structure is a flag to indicate whether the problem is rotated and added with l-p regularization term.
                    seed = min(abs(ceil(1e5*sin(1e9*sum(num_problems)))) + ...
                        abs(ceil(1e4 * sin(1e7*num_random))) + 5000 * num_random, 2^32 - 1 -num_random) + i_run;
                    rng(seed)
                    [Q,R] = qr(randn(dim));
                    Q(:, diag(R) < 0) = -Q(:, diag(R) < 0);
                    rotation_matrix = Q;
                    h = @(x) rotation_matrix * x;
                    p.objective = @(x) p.objective(h(x));
                    if isfield(parameters, "structured_factor")
                        h = @(x) parameters.structured_factor *  sum(abs(x).^ 1);
                    else
                        h = @(x) sum(abs(x).^ 1);
                    end
                    p.objective = @(x) p.objective(x) + h(x);
                end
                fval_tmp = NaN(1, num_solvers);
                if parameters.random_initial_point
                    seed = min(abs(ceil(1e5*sin(1e9*sum(num_problems)))) + ...
                        abs(ceil(1e4 * sin(1e7*num_random))) + 5000 * num_random, 2^32 - 1 -num_random) + i_run;
                    rng(seed)
                    rr = randn(size(p.x0));
                    rr = rr / norm(rr);
                    %p.x0 = p.x0 + 1 * max(1, norm(p.x0)) * rr;
                    p.x0 = p.x0 + parameters.x0_perturbation_level * max(1, norm(p.x0)) * rr;
                end
                %fprintf("%d(%d). %s\n", i_problem, i_run, p.name);
                for i_solver = 1:num_solvers
                    [fhist, fhist_perfprof] = get_fhist(p, MaxFunctionEvaluations_frec, i_solver,...
                        i_run, solvers_options, test_options);
                    if isempty(fhist)
                        fval_tmp(i_solver) = NaN;
                    else
                        fval_tmp(i_solver) = min(fhist);
                    end
                    frec(i_problem,i_solver,i_run,:) = fhist_perfprof;
                end
                fmin(i_problem, i_run) = min(fval_tmp);
            end
        end
    else
        for i_problem = 1:num_problems
            if isfield(parameters, "test_type") && strcmpi(parameters.test_type, "matcutest")
                p = macup(problem_names(1, i_problem));
            else
                problem_info = s_load(char(problem_names(i_problem)));
                p = s2mpj_wrapper(problem_info);
            end
            dim = length(p.x0);
            for i_run = 1:num_random
                % Set scaling matrix.
                if isfield(parameters, "feature") && strcmpi(parameters.feature, "badly_scaled")
                    % Badly_scaled is a flag to indicate whether the problem is badly scaled.
                    seed = min(abs(ceil(1e5*sin(1e9*sum(num_problems)))) + ...
                        abs(ceil(1e4 * sin(1e7*num_random))) + 5000 * num_random, 2^32 - 1 -num_random) + i_run;
                    rng(seed)
                    if isfield(parameters, "badly_scaled_sigma")
                        scale_matrix = diag(2.^(parameters.badly_scaled_sigma*randn(dim, 1)));
                    else
                        scale_matrix = diag(2.^(dim*randn(dim, 1)));
                    end
                    h = @(x) scale_matrix * x;
                    p.objective = @(x) p.objective(h(x));
                    p.x0 = p.x0 ./ diag(scale_matrix);
                end
                if isfield(parameters, "feature") && strcmpi(parameters.feature, "rotation_badly_scaled")
                    % Rotation_badly_scaled is a flag to indicate whether the problem is rotated and badly scaled.
                    seed = min(abs(ceil(1e5*sin(1e9*sum(num_problems)))) + ...
                        abs(ceil(1e4 * sin(1e7*num_random))) + 5000 * num_random, 2^32 - 1 -num_random) + i_run;
                    rng(seed)
                    [Q,R] = qr(randn(dim));
                    Q(:, diag(R) < 0) = -Q(:, diag(R) < 0);
                    rotation_matrix = Q;
                    h = @(x) rotation_matrix * x;
                    p.objective = @(x) p.objective(h(x));
                    p.x0 = (qr(rotation_matrix) \ eye(dim)) * p.x0;
                    scale_matrix = diag(2.^(dim*randn(dim, 1)));
                    h = @(x) scale_matrix * x;
                    p.objective = @(x) p.objective(h(x));
                    p.x0 = p.x0 ./ diag(scale_matrix);
                end
                if isfield(parameters, "feature") && (strcmpi(parameters.feature, "rotation") || ...
                        strcmpi(parameters.feature, "rotation_noisy"))
                    % Rotation is a flag to indicate whether the problem is rotated.
                    seed = min(abs(ceil(1e5*sin(1e9*sum(num_problems)))) + ...
                        abs(ceil(1e4 * sin(1e7*num_random))) + 5000 * num_random, 2^32 - 1 -num_random) + i_run;
                    rng(seed)
                    [Q,R] = qr(randn(dim));
                    Q(:, diag(R) < 0) = -Q(:, diag(R) < 0);
                    rotation_matrix = Q;
                    h = @(x) rotation_matrix * x;
                    %p.objective = @(x) p.objective(@(x) rotation_matrix * x(x));
                    p.objective = @(x) p.objective(h(x));
                    p.x0 = (qr(rotation_matrix) \ eye(dim)) * p.x0;
                end
                if isfield(parameters, "feature") && strcmpi(parameters.feature, "structured")
                    % Structured is a flag to indicate whether the problem is added with l-p regularization term.
                    if isfield(parameters, "structured_factor")
                        h = @(x) parameters.structured_factor *  sum(abs(x).^ 1);
                    else
                        h = @(x) sum(abs(x).^ 1);
                    end
                    %h = @(x) parameters.structured_factor *  sum(abs(x).^ parameters.structured_norm)^(1/parameters.structured_norm);
                    p.objective = @(x) p.objective(x) + h(x);
                end
                if isfield(parameters, "feature") && strcmpi(parameters.feature, "rotation_structured")
                    % Rotated_structure is a flag to indicate whether the problem is rotated and added with l-p regularization term.
                    seed = min(abs(ceil(1e5*sin(1e9*sum(num_problems)))) + ...
                        abs(ceil(1e4 * sin(1e7*num_random))) + 5000 * num_random, 2^32 - 1 -num_random) + i_run;
                    rng(seed)
                    [Q,R] = qr(randn(dim));
                    Q(:, diag(R) < 0) = -Q(:, diag(R) < 0);
                    rotation_matrix = Q;
                    h = @(x) rotation_matrix * x;
                    p.objective = @(x) p.objective(h(x));
                    if isfield(parameters, "structured_factor")
                        h = @(x) parameters.structured_factor *  sum(abs(x).^ 1);
                    else
                        h = @(x) sum(abs(x).^ 1);
                    end
                    p.objective = @(x) p.objective(x) + h(x);
                end
                if isfield(parameters, "rotated_badly_scaled")
                    % Rotated_badly_scaled is a flag to indicate whether the problem is rotated and badly scaled.
                    seed = min(abs(ceil(1e5*sin(1e9*sum(num_problems)))) + ...
                        abs(ceil(1e4 * sin(1e7*num_random))) + 5000 * num_random, 2^32 - 1 -num_random) + i_run;
                    rng(seed)
                    [Q,R] = qr(randn(dim));
                    Q(:, diag(R) < 0) = -Q(:, diag(R) < 0);
                    rotation_matrix = Q;
                    h = @(x) rotation_matrix * x;
                    p.objective = @(x) p.objective(h(x));
                    scale_matrix = diag(2.^(dim*randn(dim, 1)));
                    h = @(x) scale_matrix * x;
                    p.objective = @(x) p.objective(h(x));
                    p.x0 = p.x0 ./ diag(scale_matrix);
                end
                fval_tmp = NaN(1, num_solvers);
                if parameters.random_initial_point
                    seed = min(abs(ceil(1e5*sin(1e9*sum(num_problems)))) + ...
                        abs(ceil(1e4 * sin(1e7*num_random))) + 5000 * num_random, 2^32 - 1 -num_random) + i_run;
                    rng(seed)
                    rr = randn(size(p.x0));
                    rr = rr / norm(rr);
                    p.x0 = p.x0 + parameters.x0_perturbation_level * max(1, norm(p.x0)) * rr;
                end
                %fprintf("%d(%d). %s\n", i_problem, i_run, p.name);
                %fhist_tmp = cell(2, 1);
                BDS_list = ["DS", "CBDS", "PBDS", "RBDS", "PADS", "sCBDS"];
                if ~isempty(intersect(parameters.solvers_name, lower(BDS_list)))
                    % If there is a solver that we invoke existing in BDS_List, set the direction_set of input to be the same
                    % random orthogonal matrix.
                    dim = length(p.x0);
                    [direction_set_base, ~] = qr(randn(dim));
                end
                for i_solver = 1:num_solvers
                    if ismember(lower(parameters.solvers_name{i_solver}), lower(BDS_list)) ...
                            && isfield(solvers_options{i_solver}, "direction_set_type") && ...
                            strcmpi(solvers_options{i_solver}.direction_set_type, "randomized_orthogonal_matrix")
                        solvers_options{i_solver}.direction_set = direction_set_base;
                    end
                    [fhist, fhist_perfprof] = get_fhist(p, MaxFunctionEvaluations_frec, i_solver,...
                        i_run, solvers_options, test_options);
                    % if any(isnan(fhist_perfprof))
                    %     keyboard
                    % end
                    if isempty(fhist)
                        fval_tmp(i_solver) = NaN;
                    else
                        fval_tmp(i_solver) = min(fhist);
                    end
                    frec(i_problem,i_solver,i_run,:) = fhist_perfprof;
                end
                % if ~isequal(fhist_tmp{1}, fhist_tmp{2})
                %     keyboard
                % end
                fmin(i_problem, i_run) = min(fval_tmp);
            end
        end
    end

    % If parameters.fmin_type = "real-randomized", then test with the plain feature
    % should be conducted and fmin might be smaller, which makes curves
    %  of performance profile more lower.
    if strcmpi(parameters.fmin_type, "real-randomized")
        fmin_real = NaN(num_problems, 1);
        test_options.is_noisy = false;
        i_run = 1;
        if parameters.parallel == true
            parfor i_problem = 1:num_problems
                if isfield(parameters, "test_type") && strcmpi(parameters.test_type, "matcutest")
                    p = macup(problem_names(1, i_problem));
                else
                    problem_info = s_load(char(problem_names(i_problem)));
                    p = s2mpj_wrapper(problem_info);
                end
                frec_local = NaN(num_solvers, MaxFunctionEvaluations_frec);
                %fprintf("%d. %s\n", i_problem, p.name);
                for i_solver = 1:num_solvers
                    [~, fhist_perfprof] = get_fhist(p, MaxFunctionEvaluations_frec, i_solver,...
                        i_run, solvers_options, test_options);
                    frec_local(i_solver,:) = fhist_perfprof;
                end
                fmin_real(i_problem) = min(frec_local(:, :),[],"all");
            end
        else
            for i_problem = 1:num_problems
                if isfield(parameters, "test_type") && strcmpi(parameters.test_type, "matcutest")
                    p = macup(problem_names(1, i_problem));
                else
                    problem_info = s_load(char(problem_names(i_problem)));
                    p = s2mpj_wrapper(problem_info);
                end
                frec_local = NaN(num_solvers, MaxFunctionEvaluations_frec);
                %fprintf("%d. %s\n", i_problem, p.name);
                for i_solver = 1:num_solvers
                    [~, fhist_perfprof] = get_fhist(p, MaxFunctionEvaluations_frec, i_solver,...
                        i_run, solvers_options, test_options);
                    frec_local(i_solver,:) = fhist_perfprof;
                end
                fmin_real(i_problem) = min(frec_local(:, :),[],"all");
            end
        end
    end

    if strcmpi(parameters.fmin_type, "real-randomized")
        fmin_total = [fmin, fmin_real];
        fmin = min(fmin_total, [], 2);
    end

    if ~(isfield(parameters, "tuning") && parameters.tuning)
        % Plot fhist.
        compdf_location = char(fullfile(path_tests, "private", "compdf"));
        if isfield(parameters, "plot_fhist") && parameters.plot_fhist
            if num_random == 1
                outputfile = char(strcat("merged", "_", parameters.stamp_fhist, ".pdf"));
                merge_pdf(parameters.savepath, outputfile, compdf_location);
            else
                for i = 1:num_random
                    outputfile = char(strcat("merged", "_", parameters.stamp_fhist, "_", num2str(i), ".pdf"));
                    merge_pdf(parameters.savepath{i}, outputfile, compdf_location);
                end
            end
            parameters = rmfield(parameters, "savepath");
        end

        path_testdata = fullfile(path_tests, "testdata");
        path_testdata_outdir = fullfile(path_tests, "testdata", tst);

        % Make a new folder to save numerical results and source code.
        mkdir(path_testdata, tst);
        fprintf("The path of the testdata is:\n%s\n", path_testdata_outdir);
        mkdir(path_testdata_outdir, "perf");
        path_testdata_perf = fullfile(path_testdata_outdir, "perf");
        mkdir(path_testdata_perf, parameters.pdfname);
        if isfield(parameters, "log_profile") && parameters.log_profile
            log_profile = strcat(parameters.pdfname, "_", "log_perf");
            path_testdata_log_perf = fullfile(path_testdata_perf, log_profile);
            mkdir(path_testdata_log_perf);
        end
        mkdir(path_testdata_outdir, "src");
        path_testdata_src = fullfile(path_testdata_outdir, "src");
        mkdir(path_testdata_outdir, "tests");
        path_testdata_tests = fullfile(path_testdata_outdir, "tests");
        path_testdata_competitors = fullfile(path_testdata_tests, "competitors");
        mkdir(path_testdata_competitors);
        path_testdata_private = fullfile(path_testdata_tests, "private");
        mkdir(path_testdata_private);

        % Make a Txt file to store the problems that are tested and also record the dimensions of the problems that are tested.
        data_dim = zeros(1, length(problem_names));
        filePath = strcat(path_testdata_perf, "/problem_names.txt");
        fileID = fopen(filePath, 'w');
        if isfield(parameters, "test_type") && strcmpi(parameters.test_type, "matcutest")
            fprintf(fileID, '%-15s %-15s\n', 'Problem_name', 'Matcutest_dim');
        else
            fprintf(fileID, '%-15s %-15s\n', 'Problem_name', 'S2MPJ_dim');
        end
        for i_problem = 1:length(problem_names)
            if isfield(parameters, "test_type") && strcmpi(parameters.test_type, "matcutest")
                p = macup(problem_names{i_problem});
                data_dim(i_problem) = length(p.x0);
                fprintf(fileID, '%-15s %-15s\n', problem_names{i_problem}, num2str(length(p.x0)));
            else
                problem_info = s_load(char(problem_names(i_problem)));
                data_dim(i_problem) = length(problem_info.x0);
                fprintf(fileID, '%-15s %-15s\n', problem_names{i_problem}, num2str(length(problem_info.x0)));
            end
        end
        fclose(fileID);


        % Save the frequency data to a text file
        % Obtain unique values and their indices
        [uniqueValues, ~, idx] = unique(data_dim);
        % Count the number of occurrences of each unique value
        counts = histcounts(idx, 1:max(idx)+1);
        % Calculate the total count
        totalCount = sum(counts);

        % Prepare the output data
        percentage = (counts / totalCount) * 100; % Calculate the percentage
        output_data = [uniqueValues', counts', percentage']; % Combine the data

        % Only keep the rows with counts greater than 0
        output_data = output_data(counts > 0, :);

        % Sort the data by the first column
        output_data = sortrows(output_data, 1);

        % Convert the percentage to a string and add a percentage sign
        percentage_str = strcat(string(output_data(:, 3)), '%'); % Convert to string and add a percentage sign
        output_data = [output_data(:, 1:2), percentage_str]; % Merge the data with the percentage

        % Define the output path
        output_path = strcat(path_testdata_perf, "/frequency_data.txt");

        % Write the data to a text file
        writematrix(output_data, output_path, 'Delimiter', '\t', 'WriteMode', 'overwrite');

        % Display the output path
        fprintf('Frequency data saved to %s\n', output_path);

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

        % Draw performance profiles.
        % Set tolerance of convergence test in the performance profile.
        tau = parameters.tau;
        tau_length = length(tau);

        options_perf.pdfname = parameters.pdfname;
        options_perf.solvers = parameters.solvers_legend;
        options_perf.natural_stop = false;

        % Draw log-profiles if necessary.
        if isfield(parameters, "log_profile") && parameters.log_profile
            options_perf.outdir = path_testdata_log_perf;
            for l = 1:tau_length
                options_perf.tau = tau(l);
                logprof(frec, fmin, parameters.solvers_name, length(problem_names), options_perf);
            end
            outputfile = char(strcat("merged", "_", log_profile, ".pdf"));
            merge_pdf(options_perf.outdir, outputfile, compdf_location);
            movefile(fullfile(options_perf.outdir, outputfile), ...
                fullfile(path_testdata_perf, outputfile));
        end

        options_perf.outdir = fullfile(path_testdata_perf, parameters.pdfname);
        if isfield(options_perf, "tau")
            options_perf = rmfield(options_perf, "tau");
        end

        % Draw profiles.
        if parameters.is_noisy
            options_perf.feature = strcat(parameters.feature, "-", num2str(sprintf('%.1e', parameters.noise_level)));
        else
            options_perf.feature = parameters.feature;
        end
        perfdata(tau, frec, fmin, options_perf);
    end

    message = 'During the tuning process, path_testdata_perf is no needed';

    if strlength(path_testdata_perf) == 0 || all(isspace(path_testdata_perf))
        disp(message);
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