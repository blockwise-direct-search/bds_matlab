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


end
