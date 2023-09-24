function options = get_options(name_solver, options)
% GET_OPTIONS get options that needed by j-th solver on problem p.
%

bds_list = ["bds", "bds_powell"];

if any(contains(bds_list, name_solver, 'IgnoreCase', true))

    if isfield(options, "sufficient_decrease_factor_level")
        switch lower(options.sufficient_decrease_factor_level)
            case "zero"
                options.sufficient_decrease_factor = 0;
            case "negligible"
                options.sufficient_decrease_factor = 1e-16;
            case "low"
                options.sufficient_decrease_factor = 1e-8;
            case "medium"
                options.sufficient_decrease_factor = 1e-3;
            case "high"
                options.sufficient_decrease_factor = 1;
            case "excessive"
                options.sufficient_decrease_factor = 10;
            otherwise
                error("Unknown sufficient decrease factor level %s", ...
                    options.sufficient_decrease_factor_level);
        end
    end

    % Set options for matlab solvers.
elseif name_solver == "matlab_fminsearch"
    
    if isfield(options, "maxfun") && isfield(options, "StepTolerance")
        options = optimset('MaxFunEvals', options.maxfun, 'maxiter', options.maxfun,...
            'tolfun', options.StepTolerance, 'tolx', options.StepTolerance);
    elseif isfield(options, "maxfun") && ~isfield(options, "StepTolerance")    
        options = optimset('MaxFunEvals', options.maxfun, 'maxiter', options.maxfun);
    elseif ~isfield(options, "maxfun") && isfield(options, "StepTolerance")
        options = optimset('tolfun', options.StepTolerance, 'tolx', options.StepTolerance);
    end

elseif name_solver == "matlab_fminunc"
    
    if ~isfield(options, "ftarget")
        options.ftarget = -inf;
    end

    if isfield(options, "maxfun") && isfield(options, "StepTolerance")
        options = optimoptions('fminunc', 'Algorithm', 'quasi-newton', ...
            'HessUpdate', options.fminunc_type, 'MaxFunctionEvaluations',...
            options.maxfun, 'MaxIterations', options.maxfun, 'ObjectiveLimit', ...
            options.ftarget, 'StepTolerance', options.StepTolerance, ...
            'OptimalityTolerance', options.StepTolerance);
    elseif isfield(options, "maxfun") && ~isfield(options, "StepTolerance")
        options = optimoptions('fminunc', 'Algorithm', 'quasi-newton', ...
            'HessUpdate', options.fminunc_type, 'MaxFunctionEvaluations',...
            options.maxfun, 'MaxIterations', options.maxfun, 'ObjectiveLimit', ...
            options.ftarget);
    elseif ~isfield(options, "maxfun") && isfield(options, "StepTolerance")
        options = optimoptions('fminunc', 'Algorithm', 'quasi-newton', ...
            'HessUpdate', options.fminunc_type, ObjectiveLimit', options.ftarget, ...
            'StepTolerance', options.StepTolerance, 'OptimalityTolerance', ...
            options.StepTolerance);
    end

elseif name_solver == "matlab_patternsearch"
    
    if isfield(options, "maxfun") && isfield(options, "StepTolerance")
        options = optimoptions('patternsearch','MaxIterations', options.maxfun,...
            'MaxFunctionEvaluations', options.maxfun, 'FunctionTolerance',...
            options.StepTolerance, 'TolMesh', options.StepTolerance,...
            'StepTolerance', options.StepTolerance);
    elseif isfield(options, "maxfun") && ~isfield(options, "StepTolerance")
        options = optimoptions('patternsearch','MaxIterations', options.maxfun,...
            'MaxFunctionEvaluations', options.maxfun);
    elseif ~isfield(options, "maxfun") && isfield(options, "StepTolerance")
        options = optimoptions('patternsearch','MaxIterations',FunctionTolerance',...
            options.StepTolerance, 'TolMesh', options.StepTolerance,...
            'StepTolerance', options.StepTolerance);
    end

end

end
