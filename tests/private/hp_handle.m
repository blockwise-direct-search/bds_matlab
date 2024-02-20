function performance = hp_handle(value, parameters)
% perfprof_handle - Function handle for calculating the value of the
% performance profile
if ~isfield(parameters, "solvers_options")
    parameters.solvers_options = {};
end

penalty = 1000;
dist = hp_regularization(value);
x_projected = hp_projection(value);

parameters.solvers_options{1}.expand = x_projected(1);
parameters.solvers_options{1}.shrink = x_projected(2);
parameters.solvers_options{1}.reduction_factor = x_projected(3:5);

if strcmpi(parameters.solvers_name(1), "cbds")
    [frec, fmin] = hp_calculated(parameters);
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
        options_perf.tau = parameters.tau;
        performance = performance_calculated(frec, fmin, options_perf);
    end
    performance = performance + penalty * dist;
end

end


