% This script is for test.
parameters.solvers_invoke = ["cbds", "newuoa"];
parameters.problems_mindim = 1;
parameters.problems_maxdim = 1;
parameters.sufficient_decrease_factor = [0, 0];
% parameters.shuffle_period = [1, 1];
% parameters.replacement_delay = [0, 0];
%parameters.powell_factor = [0, 1e-2];
parameters.is_noisy = false;
parameters.noise_level = 1e-5;
parameters.num_random = 1;
parameters.parallel = false;
parameters.version = "now";
parameters.fmin_type = "randomized";
parameters.noise_initial_point = false;
profile_bds(parameters);