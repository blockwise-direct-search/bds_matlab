function solver_stamp = get_stamp(parameters, i)
% GET_STAMP gets the stamp of j-th solver on performance profile.

switch parameters.solvers_options{i}.solver
    case {"bds"}
        solver_stamp = upper(parameters.solvers_options{i}.Algorithm);
        if isfield(parameters.solvers_options{i}, "sufficient_decrease_factor")
            if parameters.solvers_options{i}.sufficient_decrease_factor == 0
                solver_stamp = strcat(solver_stamp, "-", ...
                    num2str(parameters.solvers_options{i}.sufficient_decrease_factor));
            elseif parameters.solvers_options{i}.sufficient_decrease_factor == eps
                solver_stamp = strcat(solver_stamp, "-", "eps");
            else
                solver_stamp = strcat(solver_stamp, "-", ...
                    int2str(int32(-log10(parameters.solvers_options{i}.sufficient_decrease_factor))));
            end
        end
    
    case {"dspd"}
        solver_stamp = "dspd";

    case {"bds_powell"}
        solver_stamp = "CBDS-Powell";

    case {"matlab_fminsearch"}
        solver_stamp = "simplex";

    case {"matlab_fminunc"}
        solver_stamp = upper(parameters.solvers_options{i}.fminunc_type);

    case {"wm_newuoa"}
        solver_stamp = "wm-newuoa";

    case {"nlopt"}
        switch parameters.solvers_options{i}.Algorithm
            case "cobyla"
                solver_stamp = "nlopt-cobyla";
            case "newuoa"
                solver_stamp = "nlopt-newuoa";
            case "bobyqa"
                solver_stamp = "nlopt-bobyqa";
        end

    case {"lam"}
        solver_stamp = "lam";
        
    case {"matlab_patternsearch"}
        solver_stamp = "patternsearch";

    case {"bfo_optimize"}
        solver_stamp = "bfo";
end

% Set solver_stamp for PRIMA family.
prima_list = ["cobyla", "uobyqa", "newuoa", "bobyqa", "lincoa"];
if any(contains(prima_list, parameters.solvers_options{i}.solver, 'IgnoreCase', true))
        solver_stamp = upper(parameters.solvers_options{i}.solver);
        if isfield(parameters.solvers_options{i}, "version")
            if strcmpi(parameters.solvers_options{i}.version, "old")
                solver_stamp = strcat(upper(parameters.solvers_options{i}), "_", "classical");
            end
        end
end

end
