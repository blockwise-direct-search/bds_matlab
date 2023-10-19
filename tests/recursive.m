function recursive(Algorithm, options)
%This function is based on https://github.com/libprima/prima/blob/main/matlab/tests/recursive.m, by Zaikun Zhang.
%RECURSIVE verifies that the solvers can be called recursively.
%

if nargin < 2
    options = struct();
end

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
orig_rng_state = rng();  % Save the current random number generator settings
rng(random_seed);

% Set the dimension of the problem
if isfield(options, 'n')
    n = options.n;
else
    n = 2;
end

% Set the recursion depth
if isfield(options, 'depth')
    depth = options.depth;
else
    depth = 3;
end

% Set up the solver
if ~isfield(options, 'compile') || options.compile
    old_directory = pwd();
    cd(fileparts(fileparts(mfilename('fullpath'))));
    setup
    cd(old_directory);
end
solver = str2func('bds');

% Define the objective function, which is based on the Rosenbrock function and recursive calls of
% the solver.
solver_options = struct();
solver_options.maxfun = min(100*n, 5e3);
solver_options.StepTolerance = 1.0e-3;
fun = @chrosen;
for i = 1 : depth
    fun = @(x) rfun(x, fun, solver, n, solver_options);
end

% Conduct the test
tic;
fprintf('\n>>>>>> Recursive test for %s starts <<<<<<\n', Algorithm);

% Call the solver
% We call the solver two times, in case something does not finish correctly during the first run.
[x, fx, exitflag, output] = solver(fun, randn(n, 1), solver_options)
[x, fx, exitflag, output] = solver(fun, randn(n, 1), solver_options)

fprintf('\n>>>>>> Recursive test for %s ends <<<<<<\n', Algorithm);
toc;

% Restore the random number generator state
rng(orig_rng_state);

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = rfun(x, fun, solver, ~, solver_options)
%RFUN defines a function of x by minimizing fun([x; y]) with respect to y in R^2 using a solver.
[~, f] = solver(@(y) fun([x; y]), randn(2, 1), solver_options);
return
