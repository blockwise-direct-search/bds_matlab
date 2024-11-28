function perf = eval_performance(solver, competitor, options)

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

    addpath(fileparts(fileparts(fileparts(mfilename('fullpath')))));
    tic;
    [~, frec, fmin] = profile(parameters);
    toc;
    perf_options = struct();
    perf_options.natural_stop = false;
    tau = options.tau;
    performances= NaN(size(tau));

    for i = 1 : length(tau)
        perf_options.tau = tau(i);
        performances(i) = performance_value(frec, fmin, perf_options);
    end
    
    perf = sum(performances.*options.weights);
    
end

