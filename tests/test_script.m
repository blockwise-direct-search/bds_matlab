clear all
parameters.test_type = "s2mpj";
parameters.parallel = true;
parameters.problem_dim = "big";
parameters.solvers_name = ["cbds", "newuoa"];
parameters.feature = "rotation_noisy_1e-3";
plot_profile(parameters);

parameters.solvers_name = ["cbds", "bfgs"];
plot_profile(parameters);

parameters.solvers_name = ["cbds", "ds"];
plot_profile(parameters);