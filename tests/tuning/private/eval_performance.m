function [perf, perf_prof_saved] = eval_performance(solver, competitor, options)

    % perf = rand();
    % return

    parameters = struct();

    parameters.tuning = true;
    parameters.parallel = false;
    parameters.problem_mindim = options.mindim;
    parameters.problem_maxdim = options.maxdim;
    parameters.num_random = options.num_random;
    parameters.test_type = options.test_type;
    parameters.feature = options.feature;
    parameters.tau = options.tau;

    parameters.solvers_name = [solver, competitor];
    parameters.solvers_options{1} = options.solver_options;
    parameters.solvers_options{1}.solver = solver;
    parameters.solvers_options{2}.solver = competitor;

    perf_prof_saved = cell(10, 1);

    addpath(fileparts(fileparts(fileparts(mfilename('fullpath')))));
    tic;
    [~, frec, fmin] = profile(parameters);
    toc;
    perf_options = struct();
    perf_options.natural_stop = false;
    performances= NaN(10, 1);

    for i = 1 : 10
        perf_options.tau = 10^(-i);
        perf_prof = performance_value(frec, fmin, perf_options);
        [performance_tuning, performance_benchmark, performances(i)] = performance_calculated(perf_prof);
        perf_prof_saved{i}.perf_prof = perf_prof;
        perf_prof_saved{i}.performance_tuning = performance_tuning;
        perf_prof_saved{i}.performance_benchmark = performance_benchmark;
        perf_prof_saved{i}.performance_diff = performances(i);
    end
    
    perf = sum(performances(1:length(options.tau))'.*options.weights);
    
end

