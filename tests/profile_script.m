parameters.test_type = "s2mpj";
parameters.parallel = true;
parameters.feature = "rotation";
parameters.problem_dim = "big";
parameters.num_random = 1;
parameters.solvers_name = ["cbds", "newuoa"];
parameters.plot_fhist = true;

plot_profile(parameters);

parameters.solvers_name = ["cbds", "bfgs"];

plot_profile(parameters);

parameters.test_type = "matcutest";

plot_profile(parameters);

parameters.solvers_name = ["cbds", "newuoa"];

plot_profile(parameters);