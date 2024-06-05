function solver_stamp = get_stamp(parameters, i)
% GET_STAMP gets the stamp of j-th solver on performance profile.
%

switch parameters.solvers_options{i}.solver
    case {"bds"}
        solver_stamp = upper(parameters.solvers_options{i}.Algorithm);
        if isfield(parameters.solvers_options{i}, "reduction_factor")
            for j = 1:length(parameters.solvers_options{i}.reduction_factor)
                if parameters.solvers_options{i}.reduction_factor(j) == 0
                    solver_stamp = strcat(solver_stamp, "-", ...
                        num2str(parameters.solvers_options{i}.reduction_factor(j)));
                elseif parameters.solvers_options{i}.reduction_factor(j) == eps
                    solver_stamp = strcat(solver_stamp, "-", "eps");
                else
                    solver_stamp = strcat(solver_stamp, "-", ...
                        int2str(int32(-log10(parameters.solvers_options{i}.reduction_factor(j)))));
                end
            end
        end

        if isfield(parameters.solvers_options{i}, "alpha_init_scaling") && parameters.solvers_options{i}.alpha_init_scaling
            solver_stamp = strcat(solver_stamp, "-", "alpha_init_scaling");
        end

        if isfield(parameters.solvers_options{i}, "forcing_function_type")
            solver_stamp = strcat(solver_stamp, "-", parameters.solvers_options{i}.forcing_function_type);
        end

        if isfield(parameters.solvers_options{i}, "forcing_function")
            if strcmp(func2str(parameters.solvers_options{i}.forcing_function), func2str(@(x)x.^2))
                solver_stamp = strcat(solver_stamp, "-", "quadratic");
            elseif strcmp(func2str(parameters.solvers_options{i}.forcing_function), func2str(@(x)x.^3))
                solver_stamp = strcat(solver_stamp, "-", "cubic");
            end
        end

        if isfield(parameters.solvers_options{i}, "forcing_function_type")
            solver_stamp = strcat(solver_stamp, "-", parameters.solvers_options{i}.forcing_function_type);
        end

        if isfield(parameters.solvers_options{i}, "shuffling_period")
            solver_stamp = strcat(solver_stamp, "-", num2str(parameters.solvers_options{i}.shuffling_period));
        end

        if isfield(parameters.solvers_options{i}, "cycling_inner")
            solver_stamp = strcat(solver_stamp, "-", num2str(parameters.solvers_options{i}.cycling_inner));
        end

        if isfield(parameters.solvers_options{i}, "replacement_delay")
            solver_stamp = strcat(solver_stamp, "-", num2str(parameters.solvers_options{i}.replacement_delay));
        end

        if isfield(parameters.solvers_options{i}, "block_indices_permuted_init")
            solver_stamp = strcat(solver_stamp, "-", "block-permuted-init");
        end

    case {"bds_norma"}
        solver_stamp = "bds_norma";

    case {"bds_previous"}
        solver_stamp = "bds_previous";

    case {"dspd"}
        solver_stamp = "dspd";
        if isfield(parameters.solvers_options{i}, "num_random_vectors")
            solver_stamp = strcat(solver_stamp, "-", num2str(parameters.solvers_options{i}.num_random_vectors));
        end

    case {"bds_powell"}
        solver_stamp = "CBDS-Powell";

    case {"fminsearch_wrapper"}
        solver_stamp = "simplex";

    case {"fminunc_wrapper"}
        solver_stamp = upper(parameters.solvers_options{i}.fminunc_type);

    case {"wm_newuoa"}
        solver_stamp = "wm-newuoa";

    case {"nlopt_wrapper"}
        switch parameters.solvers_options{i}.Algorithm
            case "cobyla"
                solver_stamp = "nlopt-cobyla";
            case "newuoa"
                solver_stamp = "nlopt-newuoa";
            case "bobyqa"
                solver_stamp = "nlopt-bobyqa";
            case "simplex"
                solver_stamp = "nlopt-simplex";
        end

    case {"imfil_wrapper"}
        solver_stamp = "imfil";

    case {"nomad_wrapper"}
        solver_stamp = "nomad";

    case {"lam"}
        solver_stamp = "lam";

    case {"patternsearch"}
        solver_stamp = "patternsearch";

    case {"bfo_wrapper"}
        solver_stamp = "bfo";

    case {"prima_wrapper"}
        solver_stamp = parameters.solvers_options{i}.Algorithm;
end

end
