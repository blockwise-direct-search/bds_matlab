function [solver_legend] = get_legend(parameters, i)
% Get the legend of solver on performance profile.

% Blockwise_direct_search
if strcmpi(parameters.solvers_invoke(i), "bds")
    solver_legend = parameters.Algorithm(i);
elseif strcmpi(parameters.solvers_invoke(i), "bds_powell")
    solver_legend = "GSDS-Powell";
end

% Matlab_fminsearch
if strcmpi(parameters.solvers_invoke(i), "matlab_fminsearch")
    solver_legend = "simplex";
end

% Matlab_fminunc
if strcmpi(parameters.solvers_invoke(i), "matlab_fminunc")
    solver_legend = parameters.fminunc_type;
end

% Prima and mnewuoa
prima_list = ["cobyla", "uobyqa", "newuoa", "bobyqa", "lincoa"];
if ~isempty(find(prima_list == parameters.solvers_invoke(i), 1))
    solver_legend = parameters.solvers_invoke(i);
    if isfield(parameters, "version")
        if strcmpi(parameters.version, "old")
            solver_legend = strcat(parameters.solvers_invoke(i), "-", "classical");
        end
    end
end

% Patternsearch
if strcmpi(parameters.solvers_invoke(i), "matlab_patternsearch")
    solver_legend = "patternsearch";
end

end
