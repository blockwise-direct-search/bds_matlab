function stress(solver, options)
% This file is cited from https://github.com/libprima/prima/blob/main/matlab/tests/stress.m, which is
% written by Zaikun ZHANG.
%
%STRESS  Stress test for the solver on problems large dimensions.
%

if nargin < 2
    options = struct();
end

% Whether to conduct a TOUGH test
tough_test = isfield(options, 'tough') && options.tough;

% Add the directory of bds.m to the path.
fullpath = mfilename('fullpath');
path_examples = fileparts(fullpath);
path_bds = fileparts(path_examples);
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
% Define the random seed by yw.
random_seed = yw;

% Set the dimension of the problem.
if isfield(options, 'n')
    n = options.n;
else
    if tough_test
        n = 1000;
    else
        n = 500;
    end
end

% Set the type of the problem.
problem_type = 'u';

% Set the type of the algorithm.
if isfield(options, "Algorithm")
    solver_options.Algorithm = options.Algorithm;
end

% Set the maxfun of the problem.
if isfield(options, "maxfun")
    solver_options.maxfun = options.maxfun;
else
    solver_options.maxfun = 500*n;
end

% The size of the array to be created. 
arraySize = [n, solver_options.maxfun];  

% Obtain the runtime instance of the current Java Virtual Machine.
rt = java.lang.Runtime.getRuntime;

% Obtain the maximum available memory size.
maxMemory = rt.maxMemory;

% Calculate the total number of bytes for the array to be created.
if isa(arraySize, 'single')
    arrayBytes = prod(arraySize) * 4;
else
    arrayBytes = prod(arraySize) * 8;
end

% Check if the array size exceeds the maximum array size limit.
if arrayBytes <= maxMemory
    solver_options.output_xhist = true;
    disp('The array size is within the range of the maximum array size limit.');
else
    disp('The array size exceeds the maximum array size limit.');
    solver_options.output_xhist = false;
end

% Set the StepTolerance of the solver.
if isfield(options, "StepTolerance")
    solver_options.StepTolerance = options.StepTolerance;
else
    solver_options.StepTolerance = 0;
end

% Generate the problem
problem = stress_problem(n, problem_type, random_seed);
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
[xval, fval, exitflag, output] = solver(problem.objective, problem.x0, solver_options)
toc;

if tough_test
    fprintf('\n>>>>>> TOUGH test ends <<<<<<\n\n');
else
    fprintf('\n>>>>>> Test ends <<<<<<\n\n');
end


return
