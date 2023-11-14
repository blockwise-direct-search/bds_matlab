function fminsearch_wrapper(FUN, x0, options)
% A wrapper of fminsearch.m.
%

% Dimension
n = numel(x0);

% Set options to an empty structure if it is not provided.
if nargin < 3
    options = struct();
end

if isfield(options, "default") && options.default
    % Set the default options.
    options = struct();
else
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

    % Set the value of StepTolerance.
    if isfield(options, "StepTolerance")
        tol = options.StepTolerance;
    else
        tol = get_default_constant("StepTolerance");
    end

    options = optimset("MaxFunEvals", maxfun, "maxiter", maxfun, "tolfun", tol, "tolx", tol);
end

% [X,FVAL,EXITFLAG,OUTPUT] = fminsearch(FUN, x0, options)
fminsearch(FUN, x0, options);

end

