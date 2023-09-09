function [parameters] = get_solvers(parameters)
% GET_SOLVERS gets the solvers that we invoke. 

solvers_num = length(parameters.solvers_options);

% RBDS - Randomized Blockwise Direct Search.
% DSPD - Direct Search based on probabilistic descent.
% CBDS - Cyclic Blockwise Direct Search.
% DS   - Direct search without blocks.
BDS_list = ["DS", "DSPD", "CBDS", "PBDS", "RBDS"];
% MATLAB_fminunc
fminunc_list = ["bfgs", "lbfgs", "dfp", "steepdesc"];

for i = 1:solvers_num 
    % If there is a solver that we invoke existing in BDS_List, set default value of Algorithm.
    if any(ismember(lower(parameters.solvers_options(i).solver), lower(BDS_list)))
        parameters.solvers_options(i).Algorithm = parameters.solvers_options(i).solver;
        parameters.solvers_options(i).solver = "bds";
    end

     % Set solver to be matlab_fminunc if it is in fminunc_list.
     if any(contains(fminunc_list, parameters.solvers_options(i).solver, 'IgnoreCase', true))
             parameters.solvers_options(i).fminunc_type = parameters.solvers_options(i).solver;
             parameters.solvers_options(i).solver = "matlab_fminunc";
     end

     % Set solver to be matlab_fminsearch if it is simplex.
     if strcmpi(parameters.solvers_options(i).solver, "simplex")
             parameters.solvers_options(i).solver = "matlab_fminsearch";
     end
end

end
