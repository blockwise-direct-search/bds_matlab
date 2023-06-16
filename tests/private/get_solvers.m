function [parameters] = get_solvers(parameters)

solvers_num = length(parameters.solvers_invoke);

% BDS(GSDS, SBDS)
BDS_list = ["GSDS", "SBDS"];
% BDS_Powell
BDS_Powell_list = "GSDS-Powell";
% Cyclic Blockwise Direct Search
CBDS_list = "CBDS";
% Direct search based on probabilistic descent (DSPD)
DSPD_list = "DSPD";
% Prima
prima_list = ["cobyla", "uobyqa", "newuoa", "bobyqa", "lincoa", "mnewuoa_wrapper"];
% Fminunc
fminunc_list = ["bfgs", "lbfgs", "dfp", "steepdesc"];

% If there is a solver called "SBDS", set default value of Algorithm.
if any(contains(parameters.solvers_invoke, "SBDS")) || any(contains(parameters.solvers_invoke, "GSDS"))...
    || any(contains(parameters.solvers_invoke, "GSDS-Powell"))
    Algorithm_list = repmat("default", 1, solvers_num);
    parameters.Algorithm = Algorithm_list;

    for i = 1:solvers_num
        if strcmp(parameters.solvers_invoke(i), "SBDS")
            parameters.Algorithm(i) = "SBDS";
        elseif strcmp(parameters.solvers_invoke(i), "GSDS")
            parameters.Algorithm(i) = "GSDS";
        end
    end
end

for i = 1:solvers_num
     % blockwise direct search
     if strcmp(parameters.solvers_invoke(i), "GSDS")...
             || strcmp(parameters.solvers_invoke(i), "SBDS" )
         parameters.solvers_invoke(i) = "bds";
     elseif ~isempty(find(BDS_Powell_list == parameters.solvers_invoke(i), 1))
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
     % CBDS
     elseif strcmp(parameters.solvers_invoke(i), "CBDS")
             parameters.solvers_invoke(i) = "bds_polling";
     % DSPD
     elseif strcmp(parameters.solvers_invoke(i), "DSPD")
                 parameters.solvers_invoke(i) = "ds_randomized";
     end
end

end
