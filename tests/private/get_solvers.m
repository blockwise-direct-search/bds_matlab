function [parameters] = get_solvers(parameters)

solvers_num = length(parameters.solvers_invoke);

% blockwise_direct_search
bds_list = ["Gauss-Seidel", "Randomized_array"];
% bds_polling
bds_polling_list = ["one", "n", "half_n", "quarter_n", "half_quarter_n"];
% ds_randomized
ds_randomized_list = ["Randomized_once", "Randomized_always"];
% Prima
prima_list = ["cobyla", "uobyqa", "newuoa", "bobyqa", "lincoa"];
% Fminunc
fminunc_list = ["bfgs", "lbfgs", "dfp", "steepdesc"];

% Set default value of parameters.blocks_strategy.
str = "none";
parameters.blocks_strategy = repmat(str, 1, solvers_num);

% If there is a solver in bds_polling_list, set default value of
% parameters.nb_generator.
if any(ismember(parameters.solvers_invoke, bds_polling_list))
    parameters.nb_generator = get_default_testparameters("nb_generator")*ones(1, solvers_num);
end

% If there is a solver in ds_randomized_list, set default value of
% parameters.randomized_strategy.
if any(ismember(parameters.solvers_invoke, ds_randomized_list))
    str = "none";
    parameters.randomized_strategy = repmat(str, 1, solvers_num);
end

for i = 1:solvers_num
     if ~isempty(find(bds_list == parameters.solvers_invoke(i), 1))
         parameters.blocks_strategy(i) = parameters.solvers_invoke(i);
         parameters.solvers_invoke(i) = "blockwise_direct_search";
     % Prima. 
     elseif ~isempty(find(prima_list == parameters.solvers_invoke(i), 1))
         
     % fminunc.    
     elseif ~isempty(find(fminunc_list == parameters.solvers_invoke(i), 1))
             parameters.fminunc_type = parameters.solvers_invoke(i);
             parameters.solvers_invoke(i) = "matlab_fminunc";
     % fminsearch        
     elseif parameters.solvers_invoke(i) == "simplex"
             parameters.solvers_invoke(i) = "matlab_fminsearch";
     % bds_polling
     elseif ~isempty(find(bds_polling_list == parameters.solvers_invoke(i), 1)) 
                 solvers_invoke = parameters.solvers_invoke(i);
                 parameters.nb_generator(i) = get_nb(solvers_invoke);
                 parameters.solvers_invoke(i) = "bds_polling";
     elseif ~isempty(find(ds_randomized_list == parameters.solvers_invoke(i), 1))             
                 parameters.randomized_strategy(i) = parameters.solvers_invoke(i);
                 parameters.solvers_invoke(i) = "ds_randomized";
     end
end

end

