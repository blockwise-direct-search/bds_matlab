function options = get_options(name_solver, options)
% GET_OPTIONS get options that needed by j-th solver on problem p.
%

bds_list = ["bds", "bds_powell"];
PRIMA_list = ["cobyla", "uobyqa", "newuoa", "bobyqa", "lincoa"];

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
elseif strcmpi(name_solver, "matlab_fminsearch")
    
    if isfield(options, "maxfun") && isfield(options, "StepTolerance")
        options = optimset('MaxFunEvals', options.maxfun, 'maxiter', options.maxfun,...
            'tolfun', options.StepTolerance, 'tolx', options.StepTolerance);
    elseif isfield(options, "maxfun") && ~isfield(options, "StepTolerance")    
        options = optimset('MaxFunEvals', options.maxfun, 'maxiter', options.maxfun);
    elseif ~isfield(options, "maxfun") && isfield(options, "StepTolerance")
        options = optimset('tolfun', options.StepTolerance, 'tolx', options.StepTolerance);
    end

elseif strcmpi(name_solver, "matlab_fminunc")
    
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

elseif strcmpi(name_solver, "matlab_patternsearch")
    
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

% Set options for PRIMA.
elseif any(contains(PRIMA_list, name_solver, 'IgnoreCase', true))

    if isfield(options, "StepTolerance")
        options.rhoend = options.StepTolerance;
    end

    if isfield(options, "alpha_init")
        options.rhobeg = options.alpha_init;
    end
    % An indicator: it can attain 0, 1, 2, 3, -1, -2, -3. Default value is
    % 0. More absolute value of iprint, more information will be printed on command
    % window. When the value of iprint is negative, no information will be
    % printed on command window and will be stored in a file.
    % options.iprint = 0;

% Set options for NLOPT.
elseif strcmpi(name_solver, "nlopt")

   if isfield(options, "StepTolerance")
        options.ftol_rel = options.StepTolerance;
        options.ftol_abs = options.StepTolerance;
   end

end

end
