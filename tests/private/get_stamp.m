function [solver_stamp] = get_stamp(parameters, j)
% GET_STAMP gets the stamp of j-th solver on performance profile.

% Set solver_stamp for BDS family.
if strcmpi(parameters.solvers_invoke(j), "bds")
    solver_stamp = upper(parameters.Algorithm(j));
    solver_stamp = strcat(solver_stamp, "_", parameters.forcing_function(j));
elseif strcmpi(parameters.solvers_invoke(j), "bds_powell")
    solver_stamp = "CBDS_Powell";
elseif strcmpi(parameters.solvers_invoke(j), "bds_cunxin")
    solver_stamp = "CBDS_Cunxin";
end

% Set solver_stamp for Matlab_fminsearch.
if strcmpi(parameters.solvers_invoke(j), "matlab_fminsearch")
    solver_stamp = "simplex";
end

% Set solver_stamp for Matlab_fminunc.
if strcmpi(parameters.solvers_invoke(j), "matlab_fminunc")
    solver_stamp = upper(parameters.fminunc_type);
end

% Set solver_stamp for PRIMA family.
prima_list = ["cobyla", "uobyqa", "newuoa", "bobyqa", "lincoa"];
if any(contains(prima_list, parameters.solvers_invoke(j), 'IgnoreCase', true))
        solver_stamp = upper(parameters.solvers_invoke(j));
        if isfield(parameters, "version")
            if strcmpi(parameters.version, "old")
                solver_stamp = strcat(upper(parameters.solvers_invoke(j)), "_", "classical");
            end
        end
end

% Set solver_stamp for patternsearch.
if strcmpi(parameters.solvers_invoke(j), "matlab_patternsearch")
    solver_stamp = "patternsearch";
end

end
