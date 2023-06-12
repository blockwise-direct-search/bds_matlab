function [solver_legend] = get_legend(parameters, i)
% Get the legend of solver on performance profile.

% Blockwise_direct_search
if strcmp(parameters.solvers_invoke(i), "blockwise_direct_search") && strcmp(parameters.blocks_strategy(i), "Gauss-Seidel")
    solver_legend = "bds-GS";
elseif strcmp(parameters.solvers_invoke(i), "blockwise_direct_search") && strcmp(parameters.blocks_strategy(i), "Randomized_array")
    solver_legend = "bds-Randomized(array)";
elseif strcmp(parameters.solvers_invoke(i), "bds_powell")
    solver_legend = "bds-GS-Powell";
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
    solver_legend = strcat("ds","_",polling_outer, num2str(parameters.cycling_outer(i)), polling_inner, ...
        num2str(parameters.cycling_inner(i)));
end

% ds_randomized
if strcmp(parameters.solvers_invoke(i), "ds_randomized")
    if parameters.randomized_strategy(i) == "Randomized_once"
        solver_legend = "ds-Randomized-once";
    elseif parameters.randomized_strategy(i) == "Randomized_always"
        solver_legend = "ds-Randomized-always";
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

% Prima and mnewuoa
prima_list = ["cobyla", "uobyqa", "newuoa", "bobyqa", "lincoa", "mnewuoa_wrapper"];
if ~isempty(find(prima_list == parameters.solvers_invoke(i), 1))
    if strcmp(parameters.solvers_invoke(i), "mnewuoa_wrapper")
        solver_legend = "mnewuoa";
    else
        solver_legend = parameters.solvers_invoke(i);
    end
end


% Patternsearch
if strcmp(parameters.solvers_invoke(i), "matlab_patternsearch")
    solver_legend = "patternsearch";
end

end

