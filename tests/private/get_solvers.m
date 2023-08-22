function [parameters] = get_solvers(parameters)
% Indicate the solvers that we invoke from its family. 

solvers_num = length(parameters.solvers_invoke);

% bds_list(RBDS, DSPD, CBDS, DS, )
% RBDS: Randomized Blockwise Direct Search
% DSPD: Direct Search based on probabilistic descent
% CBDS: Cyclic Blockwise Direct Search
% DS: Direct search without blocks
BDS_list = ["DS", "DSPD", "CBDS", "PBDS", "RBDS"];
% Prima
prima_list = ["cobyla", "uobyqa", "newuoa", "bobyqa", "lincoa"];
% Fminunc
fminunc_list = ["bfgs", "lbfgs", "dfp", "steepdesc"];

% If there is a solver called "SBDS", set default value of Algorithm.
if any(ismember(lower(parameters.solvers_invoke), lower(BDS_list)))
    Algorithm_list = repmat("default", 1, solvers_num);
    parameters.Algorithm = Algorithm_list;

    for i = 1:solvers_num
        if strcmpi(parameters.solvers_invoke(i), "ds")
            parameters.Algorithm(i) = "ds";
        elseif strcmpi(parameters.solvers_invoke(i), "dspd")
            parameters.Algorithm(i) = "dspd";
        elseif strcmpi(parameters.solvers_invoke(i), "cbds")
            parameters.Algorithm(i) = "cbds";
        elseif strcmpi(parameters.solvers_invoke(i), "pbds")
            parameters.Algorithm(i) = "pbds";
        elseif strcmpi(parameters.solvers_invoke(i), "rbds")
            parameters.Algorithm(i) = "rbds";
        end
    end
end

% If parameters.version is set to be old, then options of classical of Prima is true.
if isfield(parameters, "version")
    if strcmpi(parameters.version, "old")
        parameters.classical = true;
    end
end

for i = 1:solvers_num
     % Blockwise Direct Search
     if any(contains(BDS_list, parameters.solvers_invoke(i), 'IgnoreCase', true))
         parameters.solvers_invoke(i) = "bds";
     % Blockwise Direct Search with Powell's technique.
     elseif strcmpi(parameters.solvers_invoke(i), "Powell")
             parameters.solvers_invoke(i) = "bds_powell";
     elseif strcmpi(parameters.solvers_invoke(i), "Cunxin")
             parameters.solvers_invoke(i) = "bds_cunxin";
             parameters.Algorithm(i) = "cbds";
     % Prima.
     elseif ~isempty(find(prima_list == parameters.solvers_invoke(i), 1))

     % fminunc.
     elseif ~isempty(find(fminunc_list == parameters.solvers_invoke(i), 1))
             parameters.fminunc_type = parameters.solvers_invoke(i);
             parameters.solvers_invoke(i) = "matlab_fminunc";
     % fminsearch
     elseif parameters.solvers_invoke(i) == "simplex"
             parameters.solvers_invoke(i) = "matlab_fminsearch";
     end
end

end
