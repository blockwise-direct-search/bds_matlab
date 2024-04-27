function parameters = get_options(parameters)
% GET_OPTIONS get options that needed by solvers.
%

% Set MaxFunctionEvaluations and MaxFunctionEvaluations_factor for test_options.
if isfield(parameters, "MaxFunctionEvaluations")
    for i = 1:length(parameters.solvers_name)
        parameters.solvers_options{i}.MaxFunctionEvaluations = parameters.MaxFunctionEvaluations;
    end
end

if isfield(parameters, "MaxFunctionEvaluations_factor")
    for i = 1:length(parameters.solvers_name)
        parameters.solvers_options{i}.MaxFunctionEvaluations_factor = parameters.MaxFunctionEvaluations_factor;
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
% If the solver is fminunc_wrapper and the value of with_gradient is provided,
% we use the provided value, which implies that we are testing the case
% where fminunc_wrapper is used with and without specifying the gradient.
for i = 1:length(parameters.solvers_name)
    if strcmpi(parameters.solvers_options{i}.solver, "fminunc_wrapper")
        if ~isfield(parameters.solvers_options{i}, "with_gradient")
            parameters.solvers_options{i}.with_gradient = parameters.is_noisy;
        end        
    end
end

end
