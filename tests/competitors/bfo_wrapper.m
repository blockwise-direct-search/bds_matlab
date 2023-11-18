function bfo_wrapper(FUN, x0, options)
% A wrapper of finding the minimum of a multivariable function using derivative-free 
% method (Brute Force Optimization). For more information, please see
%https://github.com/m01marpor/BFO.
%

% Dimension
n = numel(x0);

if isfield(options, "default") && options.default
    
    %[ x, fx, msg, wrn, neval ] = bfo(FUN, x0)
    bfo(FUN, x0, 'verbosity', 'silent');

else

    if isfield(options, "StepTolerance")
        StepTolerance = options.StepTolerance;
    else
        StepTolerance = eps;
    end
    
    % Set MAXFUN to the maximum number of function evaluations.
    if isfield(options, "maxfun")
        maxeval = options.maxfun;
    else
        maxeval = get_default_constant("maxfun_dim_factor")*n;
    end
    
    %[ x, fx, msg, wrn, neval ] = bfo(FUN, x0, 'epsilon', StepTolerance, 'maxeval', maxeval)
    bfo(FUN, x0, 'epsilon', StepTolerance, 'maxeval', maxeval, 'verbosity', 'silent');
    
end

end

