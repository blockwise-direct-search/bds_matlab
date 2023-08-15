function f = tough(f, x, random_seed, noise_level, with_failure)
%This function contaminates f for the TOUGH test.
% f is the function value to be contaminated.
% x is the value of the decision variable corresponding to f; it is used when defining the random seed.
% random_seed is a seed provided by the caller in order to ensure reproducibility. The random seed
% used internally (see `rseed` below) will be defined by random_seed, f, and x.

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
elseif r > 0.6
    f = Inf;
elseif r < -0.9
    f = -1e30;
end
