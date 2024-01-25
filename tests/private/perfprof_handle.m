function performance = perfprof_handle(value, parameters, tau)
% perfprof_handle - Function handle for perfprof_calculated

if ~isfield(parameters, "solvers_options")
    parameters.solvers_options = {};
end

if value(1) < 1 || (value(2) <= 0 && value(2) >= 1) ...
        || (value(3) < 0 || value(3) > value(4) || value(3) > value(5))...
        || (value(4) <= 0 || value(4) > value(5)) || value(5) <= 0
    performance = NaN;
else
    parameters.tau = tau;
    parameters.solvers_options{1}.expand = value(1);
    parameters.solvers_options{1}.shrink = value(2);
    parameters.solvers_options{1}.reduction_factors = value(3:5);
    performance = perfprof_calculated(parameters);
end

end

