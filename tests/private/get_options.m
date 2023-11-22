function parameters = get_options(parameters)
% GET_OPTIONS get options that needed by solvers.
%

% Set maxfun and maxfun_factor for test_options.
if isfield(parameters, "maxfun")
    for i = 1:length(parameters.solvers_name)
        parameters.solvers_options{i}.maxfun = parameters.maxfun;
    end
end

if isfield(parameters, "maxfun_factor")
    for i = 1:length(parameters.solvers_name)
        parameters.solvers_options{i}.maxfun_factor = parameters.maxfun_factor;
    end
end

if isfield(parameters, "alpha_init")
    for i = 1:length(parameters.solvers_name)
        parameters.solvers_options{i}.alpha_init = parameters.alpha_init;
    end
end

if isfield(parameters, "StepTolerance")
    for i = 1:length(parameters.solvers_name)
        parameters.solvers_options{i}.StepTolerance = parameters.StepTolerance;
    end
end

if isfield(parameters, "default") && parameters.default
    for i = 1:length(parameters.solvers_name)
        parameters.solvers_options{i}.default = true;
    end
end

% Set with_gradient for fminunc. Note that when the problem is noisy, we
% should use the gradient we provide in ScalarFunction. Make sure that
% the value of with_gradient is consistent with the value in get_fhist.m.
for i = 1:length(parameters.solvers_name)
    if strcmpi(parameters.solvers_options{i}.solver, "fminunc_wrapper")
        parameters.solvers_options{i}.with_gradient = parameters.is_noisy;
    end
end

end
