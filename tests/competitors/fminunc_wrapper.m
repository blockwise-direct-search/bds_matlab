function x = fminunc_wrapper(FUN, x0, options)
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
    if isfield(options, "MaxFunctionEvaluations")
        MaxFunctionEvaluations = options.MaxFunctionEvaluations;
    else
        MaxFunctionEvaluations = get_default_constant("MaxFunctionEvaluations_dim_factor")*n;
    end

    % Set the value of StepTolerance.
    if isfield(options, "StepTolerance")
        tol = options.StepTolerance;
    else
        tol = get_default_constant("StepTolerance");
    end
    % tol = eps;

    % Set the target of the objective function.
    if isfield(options, "ftarget")
        ftarget = options.ftarget;
    else
        ftarget = get_default_constant("ftarget");
    end

    % % Set the Algorithm of fminunc.
    % if isfield(options, "fminunc_type")
    %     fminunc_type = options.fminunc_type;
    % else
    %     fminunc_type = "bfgs";
    % end

    % If and only if fminunc is invoked and the problem is noisy, with_gradient should be true, which
    % means that the gradient is provided by the user.
    with_gradient = options.with_gradient;
    
    % Set the options of fminunc. If and only if SpecifyObjectiveGradient is true,
    % fminunc will accept the gradient provided by the user. 
    % Valid value for OPTIONS parameter HessUpdate should be 'bfgs',  'lbfgs',  'dfp',  or 'steepdesc'.
    options = optimoptions("fminunc", "Algorithm", "quasi-newton", "HessUpdate", ...
        "bfgs", "MaxFunctionEvaluations", MaxFunctionEvaluations, "MaxIterations", 10^20, ...
        "ObjectiveLimit", ftarget, "StepTolerance", tol, "OptimalityTolerance", eps, ...
        'SpecifyObjectiveGradient', with_gradient);

    [x, ~, ~, ~] = fminunc(FUN, x0, options);
    %fminunc(FUN, x0, options);

end



end

