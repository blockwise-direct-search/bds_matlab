function [solver_stamp] = get_stamp(parameters, j)
% Get the stamp of solver on performance profile.

% bds
if strcmpi(parameters.solvers_invoke(j), "bds")
    solver_stamp = parameters.Algorithm(j);
elseif strcmpi(parameters.solvers_invoke(j), "bds_powell")
    solver_stamp = "CBDS_Powell";
end

% Matlab_fminsearch
if strcmpi(parameters.solvers_invoke(j), "matlab_fminsearch")
    solver_stamp = "simplex";
end

% Matlab_fminunc
if strcmpi(parameters.solvers_invoke(j), "matlab_fminunc")
    solver_stamp = parameters.fminunc_type;
end

% Prima
prima_list = ["cobyla", "uobyqa", "newuoa", "bobyqa", "lincoa"];
if ~isempty(find(prima_list == parameters.solvers_invoke(j), 1))
        solver_stamp = parameters.solvers_invoke(j);
        if isfield(parameters, "version")
            if strcmpi(parameters.version, "old")
                solver_stamp = strcat(parameters.solvers_invoke(j), "_", "classical");
            end
        end
end

% Patternsearch
if strcmpi(parameters.solvers_invoke(j), "matlab_patternsearch")
    solver_stamp = "patternsearch";
end

end
