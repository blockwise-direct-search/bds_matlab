function nomad_wrapper(fun, x0, options)
%A wrapper for NOMAD.
%

% Dimension
n = numel(x0);

% Set the default bounds.
lb = -inf(1, n);
ub = inf(1, n);

% Set MAXFUN to the maximum number of function evaluations.
if isfield(options, "MaxFunctionEvaluations_factor") && isfield(options, "MaxFunctionEvaluations")
    MaxFunctionEvaluations = min(options.MaxFunctionEvaluations_factor*n, options.MaxFunctionEvaluations);
elseif isfield(options, "MaxFunctionEvaluations_factor")
    MaxFunctionEvaluations = options.MaxFunctionEvaluations_factor*n;
elseif isfield(options, "MaxFunctionEvaluations")
    MaxFunctionEvaluations = options.MaxFunctionEvaluations;
else
    MaxFunctionEvaluations = min(get_default_constant("MaxFunctionEvaluations"), get_default_constant("MaxFunctionEvaluations_factor")*n);
end

MaxFunctionEvaluations = num2str(MaxFunctionEvaluations);

params = struct("initial_mesh_size","* 10","MAX_BB_EVAL", MaxFunctionEvaluations);

nomadOpt(fun, x0, lb, ub, params);

end

