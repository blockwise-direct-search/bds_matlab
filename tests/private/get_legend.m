function [solver_legend] = get_legend(parameters, i)
% GET_LEGEND gets the legend of i-th solver on performance profile.

% Get legend of algorithm of BDS family.
if strcmpi(parameters.solvers_invoke(i), "bds")
   solver_legend = upper(parameters.Algorithm(i));
   solver_legend = strcat(solver_legend, "-", parameters.forcing_function(i));
   %solver_legend = "our method";
elseif strcmpi(parameters.solvers_invoke(i), "bds_powell")
    solver_legend = "CBDS-Powell";
elseif strcmpi(parameters.solvers_invoke(i), "bds_cunxin")
    solver_legend = "CBDS-Cunxin";
end

% Get legend of Matlab_fminsearch.
if strcmpi(parameters.solvers_invoke(i), "matlab_fminsearch")
    solver_legend = "fminsearch";
end

% Get legend of algorithm of Matlab_fminunc family.
if strcmpi(parameters.solvers_invoke(i), "matlab_fminunc")
    solver_legend = upper(parameters.fminunc_type);
end

% Get legend of algorithm of Prima family.
prima_list = ["cobyla", "uobyqa", "newuoa", "bobyqa", "lincoa"];
if any(contains(prima_list, parameters.solvers_invoke(i), 'IgnoreCase', true))
    solver_legend = upper(parameters.solvers_invoke(i));
    if isfield(parameters, "version")
        if strcmpi(parameters.version, "old")
            solver_legend = strcat(upper(parameters.solvers_invoke(i)), "-", "classical");
        end
    end
end

% Get legend of Patternsearch.
if strcmpi(parameters.solvers_invoke(i), "wm_newuoa")
    solver_legend = "wm-newuoa";
end

% Get legend of Patternsearch.
if strcmpi(parameters.solvers_invoke(i), "matlab_patternsearch")
    solver_legend = "patternsearch";
end

end
