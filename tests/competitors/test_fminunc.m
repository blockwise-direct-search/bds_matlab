function x = test_fminunc(fun, x0, options)

    % Set options to an empty structure if it is not provided.
    if nargin < 3
        options = struct();
    end
    
    % Set MAXFUN to the maximum number of function evaluations.
    if isfield(options, "MaxFunctionEvaluations")
        MaxFunctionEvaluations = options.MaxFunctionEvaluations;
    else
        MaxFunctionEvaluations = 500 * length(x0);
    end
    
    % Set the value of StepTolerance.
    if isfield(options, "StepTolerance")
        tol = options.StepTolerance;
    else
        tol = 1e-6;
    end
    
    % Set the target of the objective function.
    if isfield(options, "ftarget")
        ftarget = options.ftarget;
    else
        ftarget = -inf;
    end
    
    % If with_gradient is specified in options.
    with_gradient = isfield(options, 'with_gradient') && options.with_gradient;
    
    % Set noise level
    if isfield(options, 'noise_level')
        noise_level = options.noise_level;
    else
        noise_level = 1e-3; % Default value if not provided
    end
    
    % Set the options of fminunc.
    options = optimoptions("fminunc", ...
        "Algorithm", "quasi-newton", ...
        "HessUpdate", "bfgs", ...
        "MaxFunctionEvaluations", MaxFunctionEvaluations, ...
        "MaxIterations", 10^20, ...
        "ObjectiveLimit", ftarget, ...
        "StepTolerance", tol, ...
        "OptimalityTolerance", eps, ...
        'SpecifyObjectiveGradient', with_gradient);
    
    % Define the inner function that computes the objective and gradient.
    fun_inner = @(x) fun_wrapper(fun, x, with_gradient, noise_level);
    
    % Call fminunc
    x = fminunc(fun_inner, x0, options);
    
    end
    
    function [f, g] = fun_wrapper(fun, x, with_gradient, noise_level)
        f = fun(x); % Call the user-defined function.
        g = []; % Initialize gradient.
        
        if with_gradient
            h = sqrt(max(abs(f), 1) * noise_level); 
            dim = length(x);
            g = NaN(dim, 1);
            V = eye(dim);
            for i = 1:dim
                f_fd = fun(x + h * V(:, i));
                g(i) = (f_fd - f) / h;
            end
    
            % Handle NaN and Inf values in the gradient.
            g(isnan(g)) = 0;
            grad_max = 10^10;
            g = min(grad_max, max(-grad_max, g)); 
        end
    end