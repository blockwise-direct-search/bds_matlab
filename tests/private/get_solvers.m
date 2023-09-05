function [parameters] = get_solvers(parameters)
% GET_SOLVERS get the solvers that we invoke. 

solvers_num = length(parameters.solvers_invoke);

% RBDS - Randomized Blockwise Direct Search.
% DSPD - Direct Search based on probabilistic descent.
% CBDS - Cyclic Blockwise Direct Search.
% DS   - Direct search without blocks.
BDS_list = ["DS", "DSPD", "CBDS", "PBDS", "RBDS"];
% Prima
prima_list = ["cobyla", "uobyqa", "newuoa", "bobyqa", "lincoa"];
% MATLAB_fminunc
fminunc_list = ["bfgs", "lbfgs", "dfp", "steepdesc"];

% If there is a solver that we invoke existing in BDS_List, set default value of Algorithm.
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
     % Set solvers_invoke to be bds if it is in BDS_list.
     if any(contains(BDS_list, parameters.solvers_invoke(i), 'IgnoreCase', true))
         parameters.solvers_invoke(i) = "bds";
     % Set solvers_invoke to be bds_powell if it is Powell.
     elseif strcmpi(parameters.solvers_invoke(i), "Powell")
             parameters.solvers_invoke(i) = "bds_powell";
    % Set solvers_invoke to be bds_cunxin if it is Cunxin.
     elseif strcmpi(parameters.solvers_invoke(i), "Cunxin")
             parameters.solvers_invoke(i) = "bds_cunxin";
             parameters.Algorithm(i) = "cbds";
     % Prima.
     elseif ~isempty(find(prima_list == parameters.solvers_invoke(i), 1))

     % Set solvers_invoke to be matlab_fminunc if it is in fminunc_list.
     elseif any(contains(fminunc_list, parameters.solvers_invoke(i), 'IgnoreCase', true))
             parameters.fminunc_type = parameters.solvers_invoke(i);
             parameters.solvers_invoke(i) = "matlab_fminunc";
    % Set solvers_invoke to be matlab_fminsearch if it is simplex.
     elseif strcmpi(parameters.solvers_invoke(i), "simplex")
             parameters.solvers_invoke(i) = "matlab_fminsearch";
     end
end

end
