parameters = struct();
parameters.expand = [1.25:0.25:2, 2.5:0.5:5];
parameters.shrink = [0.25, 0.3, 0.4, 0.5:0.05:0.85];

fprintf('Number of parameter combinations: %d\n', length(parameters.expand)*length(parameters.shrink));

solver = "cbds";
competitor = "cbds";

options = struct();
options.mindim = 6;
options.maxdim = 50;
options.test_type = "s2mpj";
options.tau = 10.^[-1:-1:-4];
options.weights = [0.3, 0.3, 0.3, 0.1];


options.feature = "plain";
fprintf('Feature:\t %s\n', options.feature);
options.num_random = 1;
plot_parameters(parameters, solver, competitor, options);   

options.feature = "noise_1e-3_no_rotation";
fprintf('Feature:\t %s\n', options.feature);
options.num_random = 3;
plot_parameters(parameters, solver, competitor, options);   

options.feature = "rotation_noisy_1e-3";
fprintf('Feature:\t %s\n', options.feature);
options.num_random = 3;
plot_parameters(parameters, solver, competitor, options);   

