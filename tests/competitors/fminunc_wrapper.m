function fminunc_wrapper(FUN, x0, options)
% A wrapper of fminsearch.m.
%

% Dimension
n = numel(x0);

% Set options to an empty structure if it is not provided.
if nargin < 3
    options = struct();
end

if isfield(options, "default") && options.default
    
    %[X,FVAL,EXITFLAG,OUTPUT] = fminunc(FUN, x0)
    fminunc(FUN, x0);

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

    % Set the target of the objective function.
    if isfield(options, "ftarget")
        ftarget = options.ftarget;
    else
        ftarget = get_default_constant("ftarget");
    end

    % Set the Algorithm of fminunc.
    if isfield(options, "fminunc_type")
        fminunc_type = options.fminunc_type;
    else
        fminunc_type = "bfgs";
    end

    options = optimoptions("fminunc", "Algorithm", "quasi-newton", "HessUpdate", ...
        fminunc_type, "MaxFunctionEvaluations", maxfun, "MaxIterations", 10^20,...
        "ObjectiveLimit", ftarget, "StepTolerance", tol, "OptimalityTolerance", eps);

    %[X,FVAL,EXITFLAG,OUTPUT] = fminunc(FUN, x0, options)
    fminunc(FUN, x0, options);

end



end

