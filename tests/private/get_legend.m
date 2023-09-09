function [solver_legend] = get_legend(parameters, i)
% GET_LEGEND gets the legend of i-th solver on performance profile.
%

switch parameters.solvers_options(i).solver
    case {"bds"}
        solver_legend = upper(parameters.solvers_options(i).Algorithm);
        solver_legend = strcat(solver_legend, "-",...
            parameters.solvers_options(i).sufficient_decrease_factor_level);
        %solver_legend = "our method";
    case {"bds_powell"}
        solver_legend = "CBDS-Powell";
    case {"matlab_fminsearch"}
        solver_legend = "fminsearch";
    case {"matlab_fminunc"}
        solver_legend = upper(parameters.solvers_options(i).fminunc_type);
    case {"wm_newuoa"}
        solver_legend = "wm-newuoa";
    case {"matlab_patternsearch"}
        solver_legend = "patternsearch";
end

% Get legend of algorithm of Prima family.
prima_list = ["cobyla", "uobyqa", "newuoa", "bobyqa", "lincoa"];
if any(contains(prima_list, parameters.solvers_options(i).solver, 'IgnoreCase', true))
    solver_legend = upper(parameters.solvers_options(i).solver);
    if isfield(parameters.solvers_options(i), "version")
        if strcmpi(parameters.solvers_options(i).version, "old")
            solver_legend = strcat(upper(parameters.solvers_options(i)), "-", "classical");
        end
    end
end



end
