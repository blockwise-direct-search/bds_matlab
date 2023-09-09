% This script is for test.
parameters.problems_mindim = 1;
parameters.problems_maxdim = 1;
parameters.is_noisy = false;
parameters.noise_level = 1e-5;
parameters.num_random = 1;
parameters.parallel = false;
parameters.random_initial_point = false;
% Establish parameters for solver 1.
struct1.solver = "bds";
struct1.Algorithm = "cbds";
struct1.sufficient_decrease_factor_level = 'low';
% Establish parameters for solver 2.
struct2.solver = "bds";
struct2.Algorithm = "pbds";
struct2.sufficient_decrease_factor_level = "high";

% Establish parameters for solver 3.
% struct3.solver = "newuoa";
% struct3.rhoend = 1e-6;
% struct3.rhobeg = 1;
% An indicator: it can attain 0, 1, 2, 3, -1, -2, -3. Default value is
% 0. More absolute value of iprint, more information will be printed on command
% window. When the value of iprint is negative, no information will be
% printed on command window and will be stored in a file.
% struct3.iprint = 0;
% struct3.output_xhist = true;

parameters.solvers_options = [struct1, struct2];

profile(parameters);