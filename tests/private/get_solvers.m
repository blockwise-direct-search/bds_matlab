function [parameters] = get_solvers(parameters)

solvers_num = length(parameters.solvers_invoke);

% BDS(GSDS, SBDS)
BDS_list = ["GSDS", "SBDS"];
% BDS_Powell
BDS_Powell_list = "GSDS-Powell";
% Cyclic Blockwise Direct Search
CBDS_list = "CBDS";
% Randomized Blockwise Direct Search
RBDS_list = "RBDS";
% Direct search based on probabilistic descent (DSPD)
DSPD_list = "DSPD";
% Prima
prima_list = ["cobyla", "uobyqa", "newuoa", "bobyqa", "lincoa", "mnewuoa_wrapper"];
% Fminunc
fminunc_list = ["bfgs", "lbfgs", "dfp", "steepdesc"];

% If there is a solver called "SBDS", set default value of Algorithm.
if any(contains(parameters.solvers_invoke, 'SBDS', 'IgnoreCase', true)) || any(contains(parameters.solvers_invoke, 'GSDS', 'IgnoreCase', true)) ...
    || any(contains(parameters.solvers_invoke, 'GSDS-Powell', 'IgnoreCase', true))
    Algorithm_list = repmat("default", 1, solvers_num);
    parameters.Algorithm = Algorithm_list;

    for i = 1:solvers_num
        if strcmpi(parameters.solvers_invoke(i), "SBDS")
            parameters.Algorithm(i) = "SBDS";
        elseif strcmpi(parameters.solvers_invoke(i), "GSDS")
            parameters.Algorithm(i) = "GSDS";
        end
    end
end

for i = 1:solvers_num
     % Blockwise Direct Search
     if strcmpi(parameters.solvers_invoke(i), "GSDS")...
             || strcmpi(parameters.solvers_invoke(i), "SBDS" )
         parameters.solvers_invoke(i) = "bds";
     % Blockwise Direct Search with Powell's technique.
     elseif strcmpi(parameters.solvers_invoke(i), "GSDS-Powell")
             parameters.solvers_invoke(i) = "bds_powell";
     elseif strcmpi(parameters.solvers_invoke(i), "RBDS")
             parameters.solvers_invoke(i) = "rbds";
     % Prima.
     elseif ~isempty(find(prima_list == parameters.solvers_invoke(i), 1))

     % fminunc.
     elseif ~isempty(find(fminunc_list == parameters.solvers_invoke(i), 1))
             parameters.fminunc_type = parameters.solvers_invoke(i);
             parameters.solvers_invoke(i) = "matlab_fminunc";
     % fminsearch
     elseif parameters.solvers_invoke(i) == "simplex"
             parameters.solvers_invoke(i) = "matlab_fminsearch";
     % CBDS
     elseif strcmpi(parameters.solvers_invoke(i), "CBDS")
             parameters.solvers_invoke(i) = "bds_polling";
     % DSPD
     elseif strcmpi(parameters.solvers_invoke(i), "DSPD")
                 parameters.solvers_invoke(i) = "ds_randomized";
     end
end

end
