function bfo_wrapper(FUN, x0, options)
% A wrapper of finding the minimum of a multivariable function using derivative-free 
% method (Brute Force Optimization). For more information, please see
%https://github.com/m01marpor/BFO.
%

% Dimension
n = numel(x0);

% Set options to an empty structure if it is not provided.
if nargin < 3
    options = struct();
end

if isfield(options, "StepTolerance")
    StepTolerance = options.StepTolerance;
else
    StepTolerance = eps;
end

% Set MAXFUN to the maximum number of function evaluations.
if isfield(options, "maxfun_factor") && isfield(options, "maxfun")
    maxeval = min(options.maxfun_factor*n, options.maxfun);
elseif isfield(options, "maxfun_factor")
    maxeval = options.maxfun_factor*n;
elseif isfield(options, "maxfun")
    maxeval = options.maxfun;
else
    maxeval = min(get_default_constant("maxfun"), get_default_constant("maxfun_factor")*n);
end

bfo(FUN, x0, 'epsilon', StepTolerance, 'maxeval', maxeval);

end

