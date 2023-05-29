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
    if parameters.polling_outer(i) == "complete"
        polling_outer = "com";
    else
        polling_outer = "opp";
    end
    if parameters.polling_inner(i) == "complete"
        polling_inner = "com";
    else
        polling_inner = "opp";
    end
    solver_legend = strcat(polling_outer, num2str(parameters.cycling_outer(i)), polling_inner, ...
        num2str(parameters.cycling_inner(i)));
end

% ds_randomized
if strcmp(parameters.solvers_invoke(i), "ds_randomized")
    if parameters.randomized_strategy(i) == "Randomized_once"
        solver_legend = "Randomized-once";
    elseif parameters.randomized_strategy(i) == "Randomized_always"
        solver_legend = "Randomized-always";
    end
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

