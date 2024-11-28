function T1_bak()

    clc
    % f = TQUARTIC(ones(3,1))
    % keyboard
    % Define a custom feature that combines "perturbed_x0", "noisy", and
    % "linearly_transformed".

    p = Problem(struct('name', 'TQUARTIC', 'fun', @(x) TQUARTIC(x), 'x0', 0.1*ones(3,1)));
    
    solvers = {@cbds_test, @ds_test};
    options.feature_name = 'custom';
    options.benchmark_id = 'cbds_eps_1000n_ds_eps_1000n_TQUARTIC_5_noisy_1e-3_no_plain';
    options.problem = p;
    options.n_runs = 1;
    options.seed = 25;
    options.solver_verbose = 2;
    options.mod_fun = @mod_fun;


    benchmark(solvers, options)
end

function f = mod_fun(x, rand_stream, problem)

    f = problem.fun(x);
    f = f + max(1, abs(f)) * 1e-3 * rand_stream.randn(1);
end

function x = cbds_test(fun, x0)

    options.Algorithm = 'cbds';
    options.StepTolerance = eps;
    options.MaxFunctionEvaluations = 1000*length(x0);
    [x, ~, ~, output] = bds(fun, x0, options);
    x
    output

end

function x = ds_test(fun, x0)

    options.Algorithm = 'ds';
    options.MaxFunctionEvaluations = 1000*length(x0);
    options.StepTolerance = eps;
    [x, ~, ~, output] = bds(fun, x0, options);
    x
    output

end