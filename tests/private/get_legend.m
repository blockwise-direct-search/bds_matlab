function solver_legend = get_legend(parameters, i)
% GET_LEGEND gets the legend of i-th solver on performance profile.
%

switch parameters.solvers_options{i}.solver
    case {"bds"}
        solver_legend = upper(parameters.solvers_options{i}.Algorithm);
        if isfield(parameters.solvers_options{i}, "sufficient_decrease_factor")
            if parameters.solvers_options{i}.sufficient_decrease_factor == 0
                solver_legend = strcat(solver_legend, "-", ...
                    num2str(parameters.solvers_options{i}.sufficient_decrease_factor));
            elseif parameters.solvers_options{i}.sufficient_decrease_factor == eps
                solver_legend = strcat(solver_legend, "-", "eps");
            else
                solver_legend = strcat(solver_legend, "-", ...
                    int2str(int32(-log10(parameters.solvers_options{i}.sufficient_decrease_factor))));
            end
        end

    case {"dspd"}
            solver_legend = "DSPD";

    case {"bds_powell"}
        solver_legend = "CBDS-Powell";

    case {"fminsearch_wrapper"}
        solver_legend = "fminsearch";

    case {"fminunc_wrapper"}
        solver_legend = upper(parameters.solvers_options{i}.fminunc_type);

    case {"wm_newuoa"}
        solver_legend = "wm-newuoa";
        
    case {"nlopt_wrapper"}
        switch parameters.solvers_options{i}.Algorithm
            case "cobyla"
                solver_legend = "nlopt-cobyla";
            case "newuoa"
                solver_legend = "nlopt-newuoa";
            case "bobyqa"
                solver_legend = "nlopt-bobyqa";
        end

    case {"patternsearch"}
        solver_legend = "patternsearch";

    case {"lam"}
        solver_legend = "lam";
        if isfield(parameters.solvers_options{i}, "linesearch_type")
            solver_legend = strcat(solver_legend, "-", ...
                parameters.solvers_options{i}.linesearch_type);
        end

    case {"bfo_wrapper"}
       solver_legend = "bfo";

    case {"prima_wrapper"}
       solver_legend = parameters.solvers_options{i}.Algorithm;
       
end

end
