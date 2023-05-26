function [solver_stamp] = get_stamp(parameters, j)
% Get the legend of solver on performance profile.

% Blockwise_direct_search
if strcmp(parameters.solvers_invoke(j), "blockwise_direct_search") && strcmp(parameters.blocks_strategy(j), "Gauss-Seidel")
    solver_stamp = "GS";
elseif strcmp(parameters.solvers_invoke(j), "blockwise_direct_search") && strcmp(parameters.blocks_strategy(j), "Randomized_array")
    solver_stamp = "Randomized_array";
end

% Matlab_fminsearch
if strcmp(parameters.solvers_invoke(j), "matlab_fminsearch")
    solver_stamp = "simplex";
end

% Matlab_fminunc
if strcmp(parameters.solvers_invoke(j), "matlab_fminunc")
    solver_stamp = "bfgs";
end

% Prima
prima_list = ["cobyla", "uobyqa", "newuoa", "bobyqa", "lincoa"];
if ~isempty(find(prima_list == parameters.solvers_invoke(j), 1))
    solver_stamp = parameters.solvers_invoke(j);
end

end

