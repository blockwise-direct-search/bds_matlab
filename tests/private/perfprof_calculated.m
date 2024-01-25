function [performance] = perfprof_calculated(parameters)
% This function is used to calculate the performance value of the solvers.

if ~isfield(parameters, "solvers_name")
    parameters.solvers_name = ["cbds", "cbds"];
end

for i = 1:length(parameters.solvers_name)
    parameters.solvers_options{i}.solver = parameters.solvers_name(i);
end

% In case no solvers are input, then throw an error.
if ~isfield(parameters, "solvers_options") || length(parameters.solvers_options) ~= 2
    error("There should be exactly two solvers when we tune the hyperparameters.")
end

% Get the parameters that the test needs.
parameters = set_profile_options(parameters);

% Tell MATLAB where to find MatCUTEst.
locate_matcutest();

% Get list of problems
s.type = "u"; % Unconstrained: 'u'
if isfield(parameters, "mindim_tuning")
    s.mindim = parameters.mindim_tuning;
else
    s.mindim = 1; % Minimum of dimension
end
if isfield(parameters, "maxdim_tuning")
    s.maxdim = parameters.maxdim_tuning;
else
    s.maxdim = 5; % Maximum of dimension
end
s.blacklist = [];
%s.blacklist = [s.blacklist, {}];

if s.mindim >= 6
    s.blacklist = [s.blacklist, { 'ARGTRIGLS', 'BROWNAL', ...
        'COATING', 'DIAMON2DLS', 'DIAMON3DLS', 'DMN15102LS', ...
        'DMN15103LS', 'DMN15332LS', 'DMN15333LS', 'DMN37142LS', ...
        'DMN37143LS', 'ERRINRSM', 'HYDC20LS', 'LRA9A', ...
        'LRCOVTYPE', 'LUKSAN12LS', 'LUKSAN14LS', 'LUKSAN17LS', 'LUKSAN21LS', ...
        'LUKSAN22LS', 'MANCINO', 'PENALTY2', 'PENALTY3', 'VARDIM',
        }];
end

if isfield(parameters, "problem_names")
    problem_names = parameters.problem_names;
else
    problem_names = secup(s);
end

fprintf("We will load %d problems\n\n", length(problem_names))

% Initialize the number of solvers.
num_solvers = length(parameters.solvers_options);
% Set MaxFunctionEvaluations_frec for performance profile.
MaxFunctionEvaluations_frec = max(get_default_profile_options("MaxFunctionEvaluations"), ...
    get_default_profile_options("MaxFunctionEvaluations_dim_factor")*parameters.problem_maxdim);

% Get the number of problems.
num_problems = length(problem_names);
% Get Number of random tests(If num_random = 1, it means no random test).
num_random = parameters.num_random;
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

% If parameters.noise_initial_point is true, then the initial point will be
% selected for each problem num_random times.
% The default value of parameters.fmin_type is set to be "randomized", then there is
% no need to test without noise, which makes the curve of the performance profile
% more higher. If parallel is true, use parfor to calculate (parallel computation),
% otherwise, use for to calculate (sequential computation).
if parameters.parallel == true
    parfor i_problem = 1:num_problems
        p = macup(problem_names(1, i_problem));
        dim = length(p.x0);
        % Set scaling matrix.
        if isfield(parameters, "feature") && strcmpi(parameters.feature, "badly_scaled")
            % Badly_scaled is a flag to indicate whether the problem is badly scaled.
            dim = length(p.x0);
            if isfield(parameters, "badly_scaled_sigma")
                scale_matrix = diag(2.^(parameters.badly_scaled_sigma*randn(dim, 1)));
            else
                scale_matrix = diag(2.^(dim*randn(dim, 1)));
            end
            % scale_matrix = hilb(dim);
            h = @(x) scale_matrix * x;
            p.objective = @(x) p.objective(h(x));
            p.x0 = inv(scale_matrix) * p.x0;
        end
        for i_run = 1:num_random
            if isfield(parameters, "plot_fhist") && parameters.plot_fhist
                fhist_plot = cell(1, num_solvers);
            end
            fval_tmp = NaN(1, num_solvers);
            if parameters.random_initial_point
                rr = randn(size(p.x0));
                rr = rr / norm(rr);
                %p.x0 = p.x0 + 1 * max(1, norm(p.x0)) * rr;
                p.x0 = p.x0 + parameters.x0_perturbation_level * max(1, norm(p.x0)) * rr;
            end
            fprintf("%d(%d). %s\n", i_problem, i_run, p.name);
            for i_solver = 1:num_solvers
                [fhist, fhist_perfprof] = get_fhist(p, MaxFunctionEvaluations_frec, i_solver,...
                    i_run, solvers_options, test_options);
                if isfield(parameters, "plot_fhist") && parameters.plot_fhist
                    fhist_plot{i_solver} = fhist;
                end
                if isempty(fhist)
                    fval_tmp(i_solver) = NaN;
                else
                    fval_tmp(i_solver) = min(fhist);
                end
                frec(i_problem,i_solver,i_run,:) = fhist_perfprof;
            end
            fmin(i_problem, i_run) = min(fval_tmp);
            if isfield(parameters, "plot_fhist") && parameters.plot_fhist
                plot_fhist(dim, fhist_plot, p.name, i_run, parameters);
            end
        end
    end
else
    for i_problem = 1:num_problems
        p = macup(problem_names(1, i_problem));
        dim = length(p.x0);
        % Set scaling matrix.
        if isfield(parameters, "feature") && strcmpi(parameters.feature, "badly_scaled")
            % Badly_scaled is a flag to indicate whether the problem is badly scaled.
            if isfield(parameters, "badly_scaled_sigma")
                scale_matrix = diag(2.^(parameters.badly_scaled_sigma*randn(dim, 1)));
            else
                scale_matrix = diag(2.^(dim*randn(dim, 1)));
            end
            %scale_matrix = hilb(dim);
            h = @(x) scale_matrix * x;
            p.objective = @(x) p.objective(h(x));
            p.x0 = inv(scale_matrix) * p.x0;
        end
        for i_run = 1:num_random
            if isfield(parameters, "plot_fhist") && parameters.plot_fhist
                fhist_plot = cell(1, num_solvers);
            end
            fval_tmp = NaN(1, num_solvers);
            if parameters.random_initial_point
                rr = randn(size(p.x0));
                rr = rr / norm(rr);
                p.x0 = p.x0 + parameters.x0_perturbation_level * max(1, norm(p.x0)) * rr;
            end
            fprintf("%d(%d). %s\n", i_problem, i_run, p.name);
            for i_solver = 1:num_solvers
                [fhist, fhist_perfprof] = get_fhist(p, MaxFunctionEvaluations_frec, i_solver,...
                    i_run, solvers_options, test_options);
                if isfield(parameters, "plot_fhist") && parameters.plot_fhist
                    fhist_plot{i_solver} = fhist;
                end
                if isempty(fhist)
                    fval_tmp(i_solver) = NaN;
                else
                    fval_tmp(i_solver) = min(fhist);
                end
                frec(i_problem,i_solver,i_run,:) = fhist_perfprof;
            end
            fmin(i_problem, i_run) = min(fval_tmp);
            if isfield(parameters, "plot_fhist") && parameters.plot_fhist
                plot_fhist(dim, fhist_plot, p.name, i_run, parameters);
            end
        end
    end
end

% If parameters.fmin_type = "real-randomized", then test without noise
% should be conducted and fmin might be smaller, which makes curves
%  of performance profile more lower.
if test_options.is_noisy && strcmpi(parameters.fmin_type, "real-randomized")
    fmin_real = NaN(num_problems, 1);
    test_options.is_noisy = false;
    i_run = 1;
    if parameters.parallel == true
        parfor i_problem = 1:num_problems
            p = macup(problem_names(1, i_problem));
            % Set scaling matrix.
            if isfield(parameters, "feature") && strcmpi(parameters.feature, "badly_scaled")
                % Badly_scaled is a flag to indicate whether the problem is badly scaled.
                dim = length(p.x0);
                scale_matrix = diag(2.^(1:dim).');
                % scale_matrix = hilb(dim);
                h = @(x) scale_matrix * x;
                p.objective = @(x) p.objective(h(x));
            end
            frec_local = NaN(num_solvers, MaxFunctionEvaluations_frec);
            if parameters.random_initial_point
                rr = randn(size(x0));
                rr = rr / norm(rr);
                p.x0 = p.x0 + parameters.x0_perturbation_level * max(1, norm(p.x0)) * rr;
            end
            fprintf("%d. %s\n", i_problem, p.name);
            for i_solver = 1:num_solvers
                frec_local(i_solver,:) = get_fhist(p, MaxFunctionEvaluations_frec,...
                    i_solver, i_run, solvers_options, test_options);
            end
            fmin_real(i_problem) = min(frec_local(:, :),[],"all");
        end
    else
        for i_problem = 1:num_problems
            p = macup(problem_names(1, i_problem));
            % Set scaling matrix.
            if isfield(parameters, "feature") && strcmpi(parameters.feature, "badly_scaled")
                % Badly_scaled is a flag to indicate whether the problem is badly scaled.
                dim = length(p.x0);
                scale_matrix = diag(2.^(1:dim).');
                %scale_matrix = hilb(dim);
                h = @(x) scale_matrix * x;
                p.objective = @(x) p.objective(h(x));
            end
            frec_local = NaN(num_solvers, MaxFunctionEvaluations_frec);
            if parameters.random_initial_point
                rr = randn(size(x0));
                rr = rr / norm(rr);
                p.x0 = p.x0 + parameters.x0_perturbation_level * max(1, norm(p.x0)) * rr;
            end
            fprintf("%d. %s\n", i_problem, p.name);
            for i_solver = 1:num_solvers
                frec_local(i_solver,:) = get_fhist(p, MaxFunctionEvaluations_frec,...
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

options_perf.natural_stop = false;
if isfield(parameters, "tau")
    options_perf.tau = parameters.tau;
else
    options_perf.tau = 1e-4;
end

performance = performance_calculated(frec, fmin, options_perf);

end

