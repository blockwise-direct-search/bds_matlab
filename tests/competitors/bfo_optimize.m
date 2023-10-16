function bfo_optimize(FUN, x0, options)
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

if isfield(options, 'maxfun')
    StepTolerance = options.StepTolerance;
else
    StepTolerance = eps;
end

if isfield(options, 'maxfun')
    maxeval = options.maxfun;
else
    maxeval = 1e3*n;
end

bfo(FUN, x0, 'epsilon', StepTolerance, 'maxeval', maxeval);

end

