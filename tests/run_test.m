% This script is for test.
parameters.solvers_invoke = ["cbds", "pbds"];
parameters.problems_mindim = 1;
parameters.problems_maxdim = 5;
parameters.sufficient_decrease_factor = [0, 0];
% parameters.shuffle_period = [1, 1];
% parameters.replacement_delay = [0, 0];
cunxin_factor = [1e-1, 1e-2];
parameters.cunxin_factor = repmat({cunxin_factor}, 1, 2);
cunxin_factor_period = 4;
parameters.cunxin_factor_period = repmat(cunxin_factor_period, 1, 2);
parameters.is_noisy = false;
parameters.noise_level = 1e-5;
parameters.num_random = 1;
parameters.parallel = true;
parameters.version = "now";
parameters.fmin_type = "randomized";
parameters.noise_initial_point = true;
profile_bds(parameters);