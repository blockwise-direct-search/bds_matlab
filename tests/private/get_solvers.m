function [parameters] = get_solvers(parameters)

solvers_num = length(parameters.solvers_invoke);

% bds
bds_list = ["Gauss-Seidel", "Randomized_array"];
% bds
bds_powell_list = "Gauss-Seidel-Powell";
% bds_polling
bds_polling_list = ["one", "n", "half_n", "quarter_n", "half_quarter_n"];
% ds_randomized
ds_randomized_list = ["Randomized_once", "Randomized_always"];
% Prima
prima_list = ["cobyla", "uobyqa", "newuoa", "bobyqa", "lincoa", "mnewuoa_wrapper"];
% Fminunc
fminunc_list = ["bfgs", "lbfgs", "dfp", "steepdesc"];

% Set default value of parameters.blocks_strategy.
if ~isfield(parameters, "blocks_strategy")
    str = "none";
    parameters.blocks_strategy = repmat(str, 1, solvers_num);
end


% If there is a solver in bds_polling_list, set default value of
% parameters.nb_generator.
if ~isfield(parameters, "nb_generator") && any(ismember(parameters.solvers_invoke, bds_polling_list))
    parameters.nb_generator = get_default_testparameters("nb_generator")*ones(1, solvers_num);
end

% If there is a solver in ds_randomized_list, set default value of
% parameters.randomized_strategy.
if ~isfield(parameters, "randomized_strategy") && any(ismember(parameters.solvers_invoke, ds_randomized_list))
    str = "none";
    parameters.randomized_strategy = repmat(str, 1, solvers_num);
end

if isfield(parameters, "powell_factor_level")
    powell_factor_level = parameters.powell_factor_level;
    parameters.powell_factor = get_powell_factor(powell_factor_level)*ones(1,solvers_num);
end
% If there is a solver in bds_powell_list, set default value of
% parameters.powell_factor.
if ~isfield(parameters, "powell_factor") && any(ismember(parameters.solvers_invoke, bds_powell_list))
     parameters.powell_factor = get_default_testparameters("powell_factor");
end

for i = 1:solvers_num
     % blockwise direct search 
     if ~isempty(find(bds_list == parameters.solvers_invoke(i), 1))
         parameters.blocks_strategy(i) = parameters.solvers_invoke(i);
         parameters.solvers_invoke(i) = "bds";
     elseif ~isempty(find(bds_powell_list == parameters.solvers_invoke(i), 1))
             parameters.solvers_invoke(i) = "bds_powell";
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

