options.maxfun = 1e4;
options.expand = 2;
options.shrink = 0.5;
options.sufficient_decrease_factor = 1e-3;
options.StepTolerance = eps;
options.ftarget = -inf;
options.polling_inner = "opportunistic";
options.polling_blocks = "Gauss-Seidel";
options.cycling_inner = 1;
options.memory = true;

fullpath = mfilename('fullpath');
[path_examples,~] = fileparts(fullpath);
[path_bds, ~, ~] = fileparts(path_examples);
path_src = fullfile(path_bds, 'src');
addpath(path_src)

[x, fval, exitflag, output] = blockwise_direct_search(@DW, [0; 0.1], options)
%[x, fval, exitflag, output] = fminsearch(@rosenb, [-1; 2], options)

rmpath(path_src)


