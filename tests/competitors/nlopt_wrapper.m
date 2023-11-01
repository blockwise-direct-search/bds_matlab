function [xopt, fmin, retcode] = nlopt_wrapper(fun, x0, options)
%A wrapper for NLopt, which is a library for nonlinear local and global optimization, 
%for functions with and without gradient information. For more information, please
%see https://github.com/stevengj/nlopt.
%

% Dimension
n = numel(x0);

if isfield(options, "Algorithm")
    switch options.Algorithm
        case "cobyla"
         opt.algorithm = NLOPT_LN_COBYLA; 
        case "newuoa"
         opt.algorithm = NLOPT_LN_NEWUOA;
        case "bobyqa"
         opt.algorithm = NLOPT_LN_BOBYQA;
    end
else
    opt.algorithm = NLOPT_LN_NEWUOA;
end

opt.min_objective = fun;

if isfield(options, "stopval")
    opt.stopval = options.stopval;
else
    opt.stopval = get_default_constant("ftarget");
end

if isfield(options, "StepTolerance")
    opt.ftol_rel = options.StepTolerance;
    opt.ftol_abs = options.StepTolerance;
else
    opt.ftol_rel = eps;
    opt.ftol_abs = eps;
end

if isfield(options, 'maxtime')
    opt.maxtime = options.maxtime;
else
    opt.maxtime = 0;
end

opt.initial_step = [2 1 1]; 

% Set MAXFUN to the maximum number of function evaluations.
if isfield(options, "maxfun_factor") && isfield(options, "maxfun")
    opt.maxeval = min(options.maxfun_factor*n, options.maxfun);
elseif isfield(options, "maxfun_factor")
    opt.maxeval = options.maxfun_factor*n;
elseif isfield(options, "maxfun")
    opt.maxeval = options.maxfun;
else
    opt.maxeval = min(get_default_constant("maxfun"), get_default_constant("maxfun_factor")*n);
end

%disp("nlopt is invoked");
[xopt, fmin, retcode] = nlopt_optimize(opt, x0');
%nlopt_optimize(opt, x0');

end