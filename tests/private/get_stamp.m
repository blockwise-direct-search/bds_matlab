function [solver_stamp] = get_stamp(parameters, j)
% Get the stamp of solver on performance profile.

% bds
if strcmp(parameters.solvers_invoke(j), "bds")
    solver_stamp = parameters.Algorithm(j);
elseif strcmp(parameters.solvers_invoke(j), "bds_powell")
    solver_stamp = "GSDS_Powell";
end

% Bds_polling
if strcmp(parameters.solvers_invoke(j), "bds_polling")
    solver_stamp = "CBDS";
end

% Ds_randomized
if strcmp(parameters.solvers_invoke(j), "ds_randomized")
    solver_stamp = "DSPD";
end

% Matlab_fminsearch
if strcmp(parameters.solvers_invoke(j), "matlab_fminsearch")
    solver_stamp = "simplex";
end

% Matlab_fminunc
if strcmp(parameters.solvers_invoke(j), "matlab_fminunc")
    solver_stamp = parameters.fminunc_type;
end

% Prima
prima_list = ["cobyla", "uobyqa", "newuoa", "bobyqa", "lincoa", "mnewuoa_wrapper"];
if ~isempty(find(prima_list == parameters.solvers_invoke(j), 1))
    if strcmp(parameters.solvers_invoke(j), "mnewuoa_wrapper")
        solver_stamp = "mnewuoa";
    else
        solver_stamp = parameters.solvers_invoke(j);
    end
end

% Patternsearch
if strcmp(parameters.solvers_invoke(j), "matlab_patternsearch")
    solver_stamp = "patternsearch";
end

end
