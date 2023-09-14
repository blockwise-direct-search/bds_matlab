function plot_profile()
% This script is for test.
parameters.problems_mindim = 1;
parameters.problems_maxdim = 5;
parameters.is_noisy = false;
parameters.noise_level = "low";
parameters.num_random = 1;
parameters.parallel = false;
parameters.random_initial_point = false;
parameters.maxfun_factor = 1e3;
% Establish parameters for solver 1.
solver1.solver = "cbds";
%solver1.classical = false;
%solver1.maxfun = 1e4;
% solver1.expand = 2;
solver1.sufficient_decrease_factor = 0;
%solver1.direction = "canonical";
%solver1.accept_simple_decrease = false;

% Establish parameters for solver 2.
solver2.solver = "dspd";
%solver2.linesearch_type = "new";
%solver2.expand = 2;
%solver2.Algorithm = "cobyla";
%solver2.maxfun = 1e4;
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

parameters.solvers_options = {};

parameters.solvers_options{1} = solver1;
parameters.solvers_options{2} = solver2;

profile(parameters);