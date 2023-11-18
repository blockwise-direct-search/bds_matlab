function nlopt_wrapper(fun, x0, options)
%A wrapper for NLopt, which is a library for nonlinear local and global optimization,
%for functions with and without gradient information. For more information, please
%see https://github.com/stevengj/nlopt.
%
% Dimension
n = numel(x0);

opt.min_objective = fun;

if isfield(options, "default") && options.default

    opt.algorithm = NLOPT_LN_NEWUOA;
    
else
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
    % Set MAXFUN to the maximum number of function evaluations.
    if isfield(options, "maxfun")
        opt.maxeval = options.maxfun;
    else
        opt.maxeval = get_default_constant("maxfun_dim_factor")*n;
    end
end

nlopt_optimize(opt, x0');

end