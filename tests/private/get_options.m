function [options] = get_options(name_solver, options)
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
                error("Unkown sufficient decrease factor level %s", ...
                    options.sufficient_decrease_factor_level);
        end
    end

    % Set options for matlab solvers.
elseif name_solver == "matlab_fminsearch"
    options = optimset('MaxFunEvals', options.maxfun, 'maxiter', options.maxfun,...
        'tolfun', options.StepTolerance, 'tolx', options.StepTolerance);

elseif name_solver == "matlab_fminunc"
    options = optimoptions('fminunc', 'Algorithm', 'quasi-newton', ...
        'HessUpdate', options.fminunc_type, 'MaxFunctionEvaluations',...
        options.maxfun, 'MaxIterations', options.maxfun, 'ObjectiveLimit', ...
        options.ftarget, 'StepTolerance', options.StepTolerance, ...
        'OptimalityTolerance', options.StepTolerance);

elseif name_solver == "matlab_patternsearch"
    options = optimoptions('patternsearch','MaxIterations', options.maxfun,...
        'MaxFunctionEvaluations', options.maxfun, 'FunctionTolerance',...
        options.StepTolerance, 'TolMesh', options.StepTolerance,...
        'StepTolerance', options.StepTolerance);

end


end
