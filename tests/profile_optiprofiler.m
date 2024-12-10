function profile_optiprofiler(options)

    clc

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Example 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % solvers = {@fminsearch_test, @fminunc_test};
    % benchmark(solvers)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Example 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % solvers = {@fminsearch_test, @fminunc_test};
    % benchmark(solvers, 'noisy')

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Example 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % solvers = {@fminsearch_test, @fminunc_test};
    % options.feature_name = 'noisy';
    % options.n_runs = 5;
    % options.problem = s_load('LIARWHD');
    % options.seed = 1;
    % benchmark(solvers, options)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Example 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ~isfield(options, 'feature_name')
        error('Please provide the feature name');
    end
    if startsWith(options.feature_name, 'noisy')
        options.noise_level = 10.^(str2double(options.feature_name(end-1:end)));
        options.feature_name = 'noisy';
    end 
    if startsWith(options.feature_name, 'rotation_noisy')
        options.noise_level = 10.^(str2double(options.feature_name(end-1:end)));
        options.feature_name = 'custom';
    end
    if startsWith(options.feature_name, 'truncated')
        options.significant_digits = str2double(options.feature_name(end));
        options.feature_name = 'truncated';
    end
    if startsWith(options.feature_name, 'random_nan')
        options.rate_nan = str2double(options.feature_name(find(options.feature_name == '_', 1, 'last') + 1:end)) / 100;
        options.feature_name = 'random_nan';
    end
    if ~isfield(options, 'labels')
        error('Please provide the labels for the solvers');
    end
    % Why we remove the truncated form feature adaptive? Fminunc do not know the noise level
    % such that it can not decide the step size.
    feature_adaptive = {'noisy', 'custom'};
    if ismember('fminunc', options.labels) && ismember(options.feature_name, feature_adaptive)
        options.labels(strcmp(options.labels, 'fminunc')) = {strcat('fminunc-adaptive-1e-', int2str(int32(-log10(options.noise_level))))};
    end
    if ~isfield(options, 'n_runs')
        if strcmpi(options.feature_name, 'plain') || strcmpi(options.feature_name, 'quantized')
            options.n_runs = 1;
        else
            options.n_runs = 3;
        end
    end
    time_str = char(datetime('now', 'Format', 'yy_MM_dd_HH_mm'));
    options.silent = false;
    options.keep_pool = true;
    options.solver_verbose = 1;
    options.p_type = 'u';
    if isfield(options, 'dim')
        if strcmpi(options.dim, 'small')
            options.mindim = 2;
            options.maxdim = 5;
        elseif strcmpi(options.dim, 'big')
            options.mindim = 6;
            options.maxdim = 50;
        end
        options = rmfield(options, 'dim');
    end
    if ~isfield(options, 'mindim')
        options.mindim = 2;
    end
    if ~isfield(options, 'maxdim')
        options.maxdim = 5;
    end
    solvers = cell(1, length(options.labels));
    for i = 1:length(options.labels)
        switch options.labels{i}
            case 'fminunc-adaptive-1e-1'
                solvers{i} = @fminunc_adaptive_1;
            case 'fminunc-adaptive-1e-2'
                solvers{i} = @fminunc_adaptive_2;
            case 'fminunc-adaptive-1e-3'
                solvers{i} = @fminunc_adaptive_3;
            case 'fminunc-adaptive-1e-4'
                solvers{i} = @fminunc_adaptive_4;
            case 'fminunc'
                solvers{i} = @fminunc_test;
            case 'fminsearch'
                solvers{i} = @fminsearch_test;
            case 'ds'
                solvers{i} = @ds_test;
            case 'pbds'
                solvers{i} = @pbds_test;
            case 'cbds'
                solvers{i} = @cbds_test;
            case 'cbds-original'
                solvers{i} = @cbds_original_test;
            case 'bfo'
                solvers{i} = @bfo_test;
            case 'newuoa'
                solvers{i} = @newuoa_test;
            case 'lam'
                solvers{i} = @lam_test;
            case 'nomad'
                solvers{i} = @nomad_test;
            otherwise
                error('Unknown solver');
        end
    end
    options.benchmark_id =[strrep(options.labels{1}, '-', '_'), '_', strrep(options.labels{2}, '-', '_'),...
        '_', num2str(options.mindim), '_', num2str(options.maxdim), '_', num2str(options.n_runs)];

    switch options.feature_name
        case 'noisy'
            options.benchmark_id = [options.benchmark_id, '_', options.feature_name, '_', int2str(int32(-log10(options.noise_level))), '_no_rotation'];
        case 'custom'
            options.benchmark_id = [options.benchmark_id, '_', 'rotation_noisy', '_', int2str(int32(-log10(options.noise_level)))];
        case 'truncated'
            options.benchmark_id = [options.benchmark_id, '_', options.feature_name, '_', int2str(options.significant_digits)];
        case 'random_nan'
            options.benchmark_id = [options.benchmark_id, '_', options.feature_name, '_', int2str(int32(options.rate_nan * 100))];
    end
    if options.run_plain
        options.benchmark_id = [options.benchmark_id, '_plain'];
    end
    options.benchmark_id = [options.benchmark_id, '_', time_str];
    options.excludelist = {'DIAMON2DLS',...
            'DIAMON2D',...
            'DIAMON3DLS',...
            'DIAMON3D',...
            'DMN15102LS',...
            'DMN15102',...
            'DMN15103LS',...
            'DMN15103',...
            'DMN15332LS',...
            'DMN15332',...
            'DMN15333LS',...
            'DMN15333',...
            'DMN37142LS',...
            'DMN37142',...
            'DMN37143LS',...
            'DMN37143',...
            'ROSSIMP3_mp',...
            'BAmL1SPLS',...
            'FBRAIN3LS',...
            'GAUSS1LS',...
            'GAUSS2LS',...
            'GAUSS3LS',...
            'HYDC20LS',...
            'HYDCAR6LS',...
            'LUKSAN11LS',...
            'LUKSAN12LS',...
            'LUKSAN13LS',...
            'LUKSAN14LS',...
            'LUKSAN17LS',...
            'LUKSAN21LS',...
            'LUKSAN22LS',...
            'METHANB8LS',...
            'METHANL8LS',...
            'SPINLS',...
            'VESUVIALS',...
            'VESUVIOLS',...
            'VESUVIOULS',...
            'YATP1CLS'};

    if strcmp(options.feature_name, 'custom')
        % We need mod_x0 to make sure that the linearly transformed problem is mathematically equivalent
        % to the original problem.
        options.mod_x0 = @mod_x0;
        % We only modify mod_fun since we are dealing with unconstrained problems.
        switch options.noise_level
            case 1e-1
                options.mod_fun = @mod_fun_1;
            case 1e-2
                options.mod_fun = @mod_fun_2;
            case 1e-3
                options.mod_fun = @mod_fun_3;
            case 1e-4
                options.mod_fun = @mod_fun_4;
            otherwise
                error('Unknown noise level');
        end
        options = rmfield(options, 'noise_level');
        options.mod_affine = @mod_affine;
    end

    benchmark(solvers, options)

end

function x0 = mod_x0(rand_stream, problem)

    [Q, R] = qr(rand_stream.randn(problem.n));
    Q(:, diag(R) < 0) = -Q(:, diag(R) < 0);
    x0 = Q * problem.x0;
end

function f = mod_fun_1(x, rand_stream, problem)

    f = problem.fun(x);
    f = f + max(1, abs(f)) * 1e-1 * rand_stream.randn(1);
end

function f = mod_fun_2(x, rand_stream, problem)

    f = problem.fun(x);
    f = f + max(1, abs(f)) * 1e-2 * rand_stream.randn(1);
end

function f = mod_fun_3(x, rand_stream, problem)

    f = problem.fun(x);
    f = f + max(1, abs(f)) * 1e-3 * rand_stream.randn(1);
end

function f = mod_fun_4(x, rand_stream, problem)

    f = problem.fun(x);
    f = f + max(1, abs(f)) * 1e-4 * rand_stream.randn(1);
end

function [A, b, inv] = mod_affine(rand_stream, problem)

    [Q, R] = qr(rand_stream.randn(problem.n));
    Q(:, diag(R) < 0) = -Q(:, diag(R) < 0);
    A = Q';
    b = zeros(problem.n, 1);
    inv = Q;
end

function x = fminsearch_test(fun, x0)

    % Dimension
    n = numel(x0);

    % Set MAXFUN to the maximum number of function evaluations.
    MaxFunctionEvaluations = 500*n;

    % Set the value of StepTolerance.
    tol = 1e-6;

    options = optimset("MaxFunEvals", MaxFunctionEvaluations, "maxiter", 10^20, "tolfun", eps, "tolx", tol);    

    x = fminsearch(fun, x0, options);
    
end

function x = fminunc_test(fun, x0)

    options = struct();
    
    % Set MAXFUN to the maximum number of function evaluations.
    if isfield(options, "MaxFunctionEvaluations")
        MaxFunctionEvaluations = options.MaxFunctionEvaluations;
    else
        MaxFunctionEvaluations = 500 * length(x0);
    end
    
    % Set the value of StepTolerance.
    if isfield(options, "StepTolerance")
        tol = options.StepTolerance;
    else
        tol = 1e-6;
    end
    
    % Set the target of the objective function.
    if isfield(options, "ftarget")
        ftarget = options.ftarget;
    else
        ftarget = -inf;
    end
    
    % Set the options of fminunc.
    options = optimoptions("fminunc", ...
        "Algorithm", "quasi-newton", ...
        "HessUpdate", "bfgs", ...
        "MaxFunctionEvaluations", MaxFunctionEvaluations, ...
        "MaxIterations", 10^20, ...
        "ObjectiveLimit", ftarget, ...
        "StepTolerance", tol, ...
        "OptimalityTolerance", eps);

    x = fminunc(fun, x0, options);

end

function x = ds_test(fun, x0)

    option.Algorithm = 'ds';
    x = bds(fun, x0, option);
end

function x = pbds_test(fun, x0)

    option.Algorithm = 'pbds';
    x = bds(fun, x0, option);
    
end

function x = cbds_test(fun, x0)

    option.Algorithm = 'cbds';
    option.expand = 1.25;
    option.shrink = 0.85;
    x = bds(fun, x0, option);
    
end

function x = cbds_original_test(fun, x0)

    option.Algorithm = 'cbds';
    x = bds(fun, x0, option);
    
end

function x = bfo_test(fun, x0)

    % Dimension
    n = numel(x0);

    StepTolerance = 1e-6;
    maxeval = 500*n;

    [x, ~, ~, ~, ~] = bfo(fun, x0, 'epsilon', StepTolerance, 'maxeval', maxeval);
    
end

function x = fminunc_adaptive_1(fun, x0)

    options.with_gradient = true;
    options.noise_level = 1e-1;
    x = test_fminunc(fun, x0, options);

end

function x = fminunc_adaptive_2(fun, x0)

    options.with_gradient = true;
    options.noise_level = 1e-2;
    x = test_fminunc(fun, x0, options);

end

function x = fminunc_adaptive_3(fun, x0)

    options.with_gradient = true;
    options.noise_level = 1e-3;
    x = test_fminunc(fun, x0, options);

end

function x = fminunc_adaptive_4(fun, x0)

    options.with_gradient = true;
    options.noise_level = 1e-4;
    x = test_fminunc(fun, x0, options);

end

function x = newuoa_test(fun, x0)

    options.maxfun = 500*length(x0);
    x = newuoa(fun, x0, options);
    
end

function x = lam_test(fun, x0)

    x = lam(fun, x0);
    
end

function x = nomad_test(fun, x0)

    % Dimension
    n = numel(x0);

    % Set the default bounds.
    lb = -inf(1, n);
    ub = inf(1, n);

    % Set MAXFUN to the maximum number of function evaluations.
    MaxFunctionEvaluations = 500*n;

    params = struct('MAX_BB_EVAL', num2str(MaxFunctionEvaluations), 'max_eval',num2str(MaxFunctionEvaluations));
    [x, ~, ~, ~, ~] = nomadOpt(fun,x0,lb,ub,params);
    
end
