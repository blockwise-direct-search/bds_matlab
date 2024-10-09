clear all
initial_value = [1.5, 0.5, 0, eps, eps];
parameters.parallel = true;
parameters.problem_mindim = 6;
parameters.problem_maxdim = 200;
parameters.test_type = "matcutest";
parameters.tuning = true;
options.MaxFunctionEvaluations = 2500;
options.verbose = true;
parameters.tuning_solver = "bds";
parameters.solvers_name = ["cbds", "cbds"];
hp_tuning(initial_value', parameters, options);

options = rmfield(options, 'verbose');
options = rmfield(options, 'MaxFunctionEvaluations');
parameters.tuning_solver = "newuoa";
options.maxfun = 2500;
hp_tuning(initial_value', parameters, options);