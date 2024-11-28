parameters = struct();
parameters.expand = [1.5, 2, 2.5, 3];
parameters.shrink = [0.25, 0.5, 0.6, 0.7, 0.8, 0.9];

solver = "cbds";
competitor = "cbds";

options = struct();
options.mindim = 2;
options.maxdim = 2;
options.test_type = "s2mpj";
options.feature = "noise_1e-3_no_rotation";
options.num_random = 1;
options.tau = [0.1, 0.01, 0.001];
options.weights = [0.5, 0.3, 0.2];

plot_parameters(parameters, solver, competitor, options);   
