function prima_wrapper(fun, x0, options)
%A wrapper for PRIMA.
%

% Dimension
n = numel(x0);

if isfield(options, "StepTolerance")
    options.rhoend = options.StepTolerance;
else
    options.rhoend = get_default_constant("StepTolerance");
end

if isfield(options, "alpha_init")
    options.rhobeg = options.alpha_init;
else
    options.rhobeg = get_default_constant("alpha_init");
end
% An indicator: it can attain 0, 1, 2, 3, -1, -2, -3. Default value is
% 0. More absolute value of iprint, more information will be printed on command
% window. When the value of iprint is negative, no information will be
% printed on command window and will be stored in a file.
% options.iprint = 0;

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

options.maxfun = maxfun;
if isfield(options, "maxfun_factor")
    options = rmfield(options, "maxfun_factor");
end

if isfield(options, "Algorithm")
    Algorithm = options.Algorithm;
else
    Algorithm = "newuoa";
end

solver = str2func(Algorithm);

% Turn off warning to save computation resource.
warnoff(Algorithm);

solver(fun, x0, options);

% Turn off warning is a very dangerous thing. So it must be set a loop to
% turn on after ending the computation.
warnoff(Algorithm);

end

