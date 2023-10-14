function parallel(Algorithm, options)
%This file is cited from https://github.com/libprima/prima/blob/main/matlab/tests/parallel.m, which
%is written by Zaikun Zhang.
%PARALLEL verifies that the solvers can be called in parallel.

% Turn off the warning about the debug mode.
orig_warning_state = warning;
warning('off', ['bds', ':Debug']);

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
    n = 20;
end

% Set the number of parallel runs
if isfield(options, 'np')
    np = options.np;
else
    np = 32;
end

% Set up the solver
if ~isfield(options, 'compile') || options.compile
    old_directory = pwd();
    cd(fileparts(fileparts(mfilename('fullpath'))));
    setup bds
    cd(old_directory);
end
solver = str2func('bds');

% Conduct the test
tic;
fprintf('\n>>>>>> Parallel test for %s starts <<<<<<\n', Algorithm);

% Call the solver
solver_options = struct();
solver_options.Algorithm = Algorithm;
solver_options.StepTolerance = 1.0e-8;

% We conduct two parallel tests, in case something does not finish correctly during the first run.
for i = 1 : 2

    ticBytes(gcp)

    parfor i = 1:np
        fprintf('\n>>>>>> Parallel test for %s, %d-th run <<<<<<\n', Algorithm, i);
        test(solver, n, solver_options, random_seed + i);
    end

    tocBytes(gcp)

end

fprintf('\n>>>>>> Parallel test for %s ends <<<<<<\n', Algorithm);
toc;

% Restore the random number generator state
rng(orig_rng_state);
% Restore the warning state
warning(orig_warning_state);

return


function test(solver, n, solver_options, random_seed)
rng(random_seed);
shift = randn(n, 1);
fun = @(x) chrosen(x + shift);
x0 = randn(n, 1);
[~, g0] = chrosen(x0 + shift);
[x, fx, exitflag, output] = solver(fun, x0, solver_options)
[~, gx] = chrosen(x + shift);

whos

norm(gx) / norm(g0)
if strcmpi(solver_options.Algorithm, 'ds')
    assert(norm(gx) < 1.0e-2 * norm(g0), 'X is close to stationary');
else
    assert(norm(gx) < 1.0e-3 * norm(g0), 'X is close to stationary');
end

return
