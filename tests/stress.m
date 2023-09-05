function [] = stress(solver, options)
%STRESS  Stress test for the solver on problems large dimensions.
%

if nargin < 2
    options = struct();
end
    
    % Whether to conduct a TOUGH test
    tough_test = isfield(options, 'tough') && options.tough;
    
    % Add the directory of bds.m to the path.
    fullpath = mfilename('fullpath');
    [path_examples,~] = fileparts(fullpath);
    [path_bds, ~, ~] = fileparts(path_examples);
    path_src = fullfile(path_bds, 'src');
    addpath(path_src);
    
    % Set up the solver.
    solver = str2func(solver);
    
    % Set the random seed. We ALTER THE SEED WEEKLY to test the solvers as much as possible.
    if isfield(options, 'yw')
        yw = options.yw;
    elseif isfield(options, 'seed')
        yw = options.seed;
    else
        yw = year_week('Asia/Shanghai');
    end
    fprintf('\nYW = %d\n', yw);
    % Define the random seed by yw
    random_seed = yw;
    
    % Set the dimension of the problem.
    if isfield(options, 'n')
        n = options.n;
    else
        if tough_test
            n = 200;
        else
            n = 500;
        end
    end
    
    % Set the type of the problem.
    problem_type = 'u';
    
    % Set the options for the test.
    test_options = struct();
    test_options.maxfun = 1000 * n;
    test_options.alpha_init = 1;
    test_options.StepTolerance = eps;
    if ~isfield(options, "Algorithm")
        options.Algorithm = "cbds";
    end
    test_options.Algorithm = options.Algorithm;
    fprintf('\n>>>>>> test_options =');
    test_options
    
    % Generate the problem
    problem = stress_problem(n, problem_type, random_seed);
    problem.options = test_options;
    original_problem = problem;
    if tough_test
        problem = tough(original_problem, random_seed);
    end
    
    % Conduct the test.
    if tough_test
        fprintf('\n>>>>>> TOUGH test starts <<<<<<\n');
    else
        fprintf('\n>>>>>> Test starts <<<<<<\n');
    end
    
    tic;
    exception = [];
    try
        solver(problem.objective, problem.x0, problem.options);
    catch exception
    end
    toc;    
    
    if ~isempty(exception)
        rethrow(exception);
    end
    
    if tough_test
        fprintf('\n>>>>>> TOUGH test ends <<<<<<\n\n');
    else
        fprintf('\n>>>>>> Test ends <<<<<<\n\n');
    end
    
    
    return