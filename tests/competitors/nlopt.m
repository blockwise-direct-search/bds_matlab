function nlopt(fun, x0, options)

% Dimension
n = numel(x0);

opt.algorithm = NLOPT_LN_COBYLA;
%opt.algorithm = options.Algorithm;
opt.min_objective = @(x)objective(fun,x);

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

if isfield(options, 'maxfun')
    opt.maxeval = options.maxfun;
else
    opt.maxeval = 1e3*n;
end

nlopt_optimize(opt, x0');

end
