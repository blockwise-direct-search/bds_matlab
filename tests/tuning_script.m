clear all
initial_value = [1.25, 0.8];
parameters.parallel = true;
parameters.problem_mindim = 1;
parameters.problem_maxdim = 5;
parameters.test_type = "s2mpj";
parameters.tuning = true;
parameters.feature = "noise_1e-3_no_rotation";
% parameters.feature = "plain";
parameters.tuning_solver = "bobyqa";
parameters.min_precision = 4;
parameters.solvers_name = ["cbds", "cbds"];
options.output_xhist = true;
options.iprint = 3;
% options.MaxFunctionEvaluations = 1;
hp_tuning(initial_value', parameters, options);