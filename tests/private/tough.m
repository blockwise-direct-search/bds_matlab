function tough_problem = tough(problem, random_seed, noise_level, with_failure)
    % This file is cited from https://github.com/libprima/prima/blob/main/matlab/tests/private/tough.m, which is
    % written by Zaikun Zhang.
    %This function prepares the TOUGH version of a given problem.
    % problem: a structure defining the original problem
    % random_seed: a seed provided by the caller in order to ensure reproducibility
    % noise_level: level of the noise
    % with_failure: whether to fail the objective and constraint evaluation randomly
    
    if nargin < 3
        noise_level = 2.0e-1;  % The noise level.
    end
    if nargin < 4
        with_failure = true;  % Whether to fail the function evaluation randomly.
    end
    
    % Set the random seed
    orig_rng_state = rng();
    rng(random_seed);
    
    % Copy the problem options
    if isfield(problem, 'options')
        tough_problem.options = problem.options;
    else
        tough_problem.options = [];
    end
    
    % Set the starting point
    x0 = problem.x0;
    n = length(x0);
    tough_problem.x0 = x0 + noise_level * max(1, abs(x0)) .* randn(n,1);
    
    % Set the objective function
    tough_problem.objective = @(x) tough_feval(problem.objective, x, random_seed, noise_level, with_failure);
    
    % Restore the random seed
    rng(orig_rng_state);
    
    % `tough` ends here
    return
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function f = tough_feval(func, x, random_seed, noise_level, with_failure)
    %This function evaluates the function func at x for the TOUGH test.
    if nargin < 4
        noise_level = 2e-1;
    end
    if nargin < 5
        with_failure = true;
    end
    f = func(x);
    f = contaminate(f, x, random_seed, noise_level, with_failure);
    
    % `tough_feval` ends here
    return
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function f = contaminate(f, x, random_seed, noise_level, with_failure)
    %This function contaminates f for the TOUGH test.
    % f is the function value to be contaminated.
    % x is the value of the decision variable corresponding to f.
    % The random seed used internally (see `rseed` below) will be defined by random_seed, f, and x.
    
    if nargin < 4
        noise_level = 2e-1;  % The noise level.
    end
    if nargin < 5
        with_failure = true;  % Whether to fail the function evaluation randomly.
    end
    
    % Set the random seed.
    orig_rng_state = rng();
    rseed = max(0, min(2^32 - 1, random_seed + sum(num2str(f, 16)) + sum(num2str(x, 16), 'all')));
    rng(rseed);
    
    % Contaminate f. The value will be further modified below.
    f = f * (1 + noise_level * randn);
    
    % Generate a random number to decide how to modify f.
    r = 2 * rand - 1;
    
    % Restore the random seed. Do this before the possible invocation of `error`.
    rng(orig_rng_state);
    
    % Modify the value of f to make it "tough".
    if r > 0.9
        if with_failure
            error('Function evaluation fails!');
        else
            f = NaN;
        end
    elseif r > 0.8
        f = NaN;
    elseif r > 0.7
        f = Inf;
    elseif r < -0.9
        f = -1e30;
    end
    
    % `contaminate` ends here
    return