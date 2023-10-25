function nomad_wrapper(fun, x0, options)
%A wrapper for NOMAD.
%

% Dimension
n = numel(x0);

% Set the default bounds.
lb = -inf(1, n);
ub = inf(1, n);

% Set MAXFUN to the maximum number of function evaluations.
if isfield(options, "maxfun_factor") && isfield(options, "maxfun")
    maxfun = min(options.maxfun_factor*n, options.maxfun);
elseif isfield(options, "maxfun_factor")
    maxfun = options.maxfun_factor*n;
elseif isfield(options, "maxfun")
    maxfun = options.maxfun;
else
    maxfun = min(get_default_constant("maxfun"), get_default_constant("maxfun_factor")*n);
end

maxfun = num2str(maxfun);

params = struct("initial_mesh_size","* 10","MAX_BB_EVAL", maxfun);

nomadOpt(fun, x0, lb, ub, params);

end

