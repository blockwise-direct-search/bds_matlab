function [solver_legend] = get_legend(parameters, i)
% Get the legend of solver on performance profile.

% Blockwise_direct_search
if strcmp(parameters.solvers_invoke(i), "blockwise_direct_search") && strcmp(parameters.blocks_strategy(i), "Gauss-Seidel")
    solver_legend = "GS";
elseif strcmp(parameters.solvers_invoke(i), "blockwise_direct_search") && strcmp(parameters.blocks_strategy(i), "Randomized_array")
    solver_legend = "Randomized(array)";
end

% Bds_polling
if strcmp(parameters.solvers_invoke(i), "bds_polling")
    solver_legend = "bds(polling)";
end

% Matlab_fminsearch
if strcmp(parameters.solvers_invoke(i), "matlab_fminsearch")
    solver_legend = "simplex";
end

% Matlab_fminunc
if strcmp(parameters.solvers_invoke(i), "matlab_fminunc")
    solver_legend = parameters.fminunc_type;
end

% Prima
prima_list = ["cobyla", "uobyqa", "newuoa", "bobyqa", "lincoa"];
if ~isempty(find(prima_list == parameters.solvers_invoke(i), 1))
    solver_legend = parameters.solvers_invoke(i);
end

% Patternsearch
if strcmp(parameters.solvers_invoke(i), "matlab_patternsearch")
    solver_legend = "patternsearch";
end

end

