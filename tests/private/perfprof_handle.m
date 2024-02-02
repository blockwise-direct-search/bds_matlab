function performance = perfprof_handle(value, parameters)
% perfprof_handle - Function handle for calculating the value of the
% performance profile

if ~isfield(parameters, "solvers_options")
    parameters.solvers_options = {};
end

if strcmpi(parameters.solvers_name(1), "cbds")
    if value(1) < 1 || (value(2) <= 0 && value(2) >= 1) ...
            || (value(3) < 0 || value(3) > value(4) || value(3) > value(5))...
            || (value(4) <= 0 || value(4) > value(5)) || value(5) <= 0
        performance = NaN;
    else
        parameters.solvers_options{1}.expand = value(1);
        parameters.solvers_options{1}.shrink = value(2);
        parameters.solvers_options{1}.reduction_factor = value(3:5);
        [frec, fmin] = perfprof_calculated(parameters);
        options_perf.natural_stop = false;
        if length(parameters.tau) > 1
            num_tau = length(parameters.tau);
            multi_performance = NaN(num_tau, 1);
            for i = 1 : num_tau
                options_perf.tau = parameters.tau(i);
                multi_performance(i) = performance_calculated(frec, fmin, options_perf);
            end
            performance = max(multi_performance);
        else
            performance = performance_calculated(frec, fmin, options_perf);
        end


    end
end
end


