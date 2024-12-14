function x = nomad_wrapper(fun, x0, options)
%A wrapper for NOMAD.
%

% Dimension:
n = numel(x0);

% Set the default bounds.
lb = -inf(n, 1);
ub = inf(n, 1);

% Set MAXFUN to the maximum number of function evaluations.
if isfield(options, "MaxFunctionEvaluations")
    MaxFunctionEvaluations = options.MaxFunctionEvaluations;
else
    MaxFunctionEvaluations = get_default_constant("MaxFunctionEvaluations_dim_factor")*n;
end

params = struct('MAX_BB_EVAL', num2str(MaxFunctionEvaluations), 'max_eval',num2str(MaxFunctionEvaluations));

% As of NOMAD version 4.4.0 and OptiProfiler commit 24d8cc0, the following line is 
% necessary. Otherwise, NOMAD will throw an error, complaining that the blackbox 
% evaluation fails. This seems to be because OptiProfiler wraps the function 
% handle in a way that NOMAD does not expect: NOMAD expects a function handle 
% `fun` with the signature fun(x), where x is a column vector, while OptiProfiler 
% produces one with the signature @(varargin)featured_problem.fun(varargin{:}).
fun = @(x) fun(x(:));

[x, ~, ~, ~, ~] = nomadOpt(fun,x0,lb,ub,params);

end