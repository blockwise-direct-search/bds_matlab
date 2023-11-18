function fminsearch_wrapper(FUN, x0, options)
% A wrapper of fminsearch.m.
%

% Dimension
n = numel(x0);

if isfield(options, "default") && options.default
    % Set the default options.
    fminsearch(FUN, x0);
else
    % Set MAXFUN to the maximum number of function evaluations.
    if isfield(options, "maxfun")
        maxfun = options.maxfun;
    else
        maxfun = get_default_constant("maxfun_dim_factor")*n;
    end

    % Set the value of StepTolerance.
    if isfield(options, "StepTolerance")
        tol = options.StepTolerance;
    else
        tol = get_default_constant("StepTolerance");
    end

    options = optimset("MaxFunEvals", maxfun, "maxiter", 10^20, "tolfun", eps, "tolx", tol);

    % [X,FVAL,EXITFLAG,OUTPUT] = fminsearch(FUN, x0, options)
    fminsearch(FUN, x0, options);
    
end



end

