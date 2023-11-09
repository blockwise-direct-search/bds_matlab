function plot_profile(parameters)
% This script is for the performance profile of the solvers on some CUTEst problems. The basic input of parameters
% should be parameters.solvers_name = ["solver1", "solver2" ...]. Additionally, if the user wants to add
% options for the solvers, the input should be parameters.solvers_options = {options1, options2 ...},
% where options are structures.  
% 

if ~isfield(parameters, "solvers_options")
    parameters.solvers_options = {};
end

for i = 1:length(parameters.solvers_name)
    parameters.solvers_options{i}.solver = parameters.solvers_name(i);
end

profile(parameters);