function solver_stamp = get_stamp(parameters, i)
% GET_STAMP gets the stamp of j-th solver on performance profile.

switch parameters.solvers_options{i}.solver
    case {"bds"}
        solver_stamp = upper(parameters.solvers_options{i}.solver);
        if isfield(parameters.solvers_options{i}, "sufficient_decrease_factor_level")
            solver_stamp = strcat(solver_stamp, "_",...
                parameters.solvers_options{i}.sufficient_decrease_factor_level);
        end
        %solver_legend = "our method";
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
                solver_stamp = "nlopt_cobyla";
            case "newuoa"
                solver_stamp = "nlopt_newuoa";
            case "bobyqa"
                solver_stamp = "nlopt_bobyqa";
        end
    case {"lam"}
        solver_stamp = "lam";
    case {"matlab_patternsearch"}
        solver_stamp = "patternsearch";
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
