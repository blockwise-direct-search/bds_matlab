function nlopt(fun, x0, options)

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
    opt.stopval = -inf;
end

if isfield(options, 'ftol_rel')
    opt.ftol_rel = options.ftol_rel;
else
    opt.ftol_rel = eps;
end

if isfield(options, 'ftol_abs')
    opt.ftol_abs = options.ftol_abs;
else
    opt.ftol_abs = eps;
end

if isfield(options, 'maxtime')
    opt.maxtime = options.maxtime;
else
    if n <= 5
        opt.maxtime = 20;
    else 
        opt.maxtime = 30;
    end
end

if isfield(options, 'maxfun')
    opt.maxeval = options.maxfun;
else
    opt.maxeval = 1e3*n;
end

nlopt_optimize(opt, x0');

end
