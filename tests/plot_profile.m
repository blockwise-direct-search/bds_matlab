function plot_profile(parameters)
% This script is for test.
parameters.problems_dim = "small";
parameters.maxfun_factor = 1e3;
parameters.alpha_init = 1;
parameters.StepTolerance = eps;

if ~isfield(parameters, "solvers_options")
    parameters.solvers_options = {};
end

for i = 1:length(parameters.solvers_name)
    parameters.solvers_options{i}.solver = parameters.solvers_name(i);
end

profile(parameters);