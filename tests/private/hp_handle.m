function performance = hp_handle(value, parameters)
% perfprof_handle - Function handle for calculating the value of the
% performance profile
if ~isfield(parameters, "solvers_options")
    parameters.solvers_options = {};
end

parameters.solvers_options{1}.expand = value(1);
parameters.solvers_options{1}.shrink = value(2);
parameters.solvers_options{1}.reduction_factor = value(3:5);

if strcmpi(parameters.solvers_name(1), "cbds")
    if value(1) < 1 || (value(2) <= 0 && value(2) >= 1) ...
            || (value(3) < 0 || value(3) > value(4) || value(3) > value(5))...
            || (value(4) <= 0 || value(4) > value(5)) || value(5) <= 0
        performance = NaN;
    else
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
    end
end

if strcmpi(parameters.solvers_name(1), "pbds")
    if value(1) < 1 || (value(2) <= 0 && value(2) >= 1) ...
            || (value(3) < 0 || value(3) > value(4) || value(3) > value(5))...
            || (value(4) <= 0 || value(4) > value(5)) || value(5) <= 0 ...
            || value(6) < 1
        performance = NaN;
    else
        options_perf.natural_stop = false;
        if isinteger(value(6))
            parameters.solvers_options{1}.permuting_period = value(6);
            [frec, fmin] = hp_calculated(parameters);
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
        else
            truncated_performance = NaN(2, 1);
            for i = 1 : 2
                parameters.solvers_options{1}.permuting_period = floor(value(6) - 1) + i;
                [frec, fmin] = hp_calculated(parameters);
                if length(parameters.tau) > 1
                    num_tau = length(parameters.tau);
                    multi_performance = NaN(num_tau, 1);
                    for j = 1 : num_tau
                        options_perf.tau = parameters.tau(i);
                        multi_performance(j) = performance_calculated(frec, fmin, options_perf);
                    end
                    truncated_performance(i) = max(multi_performance);
                else
                    options_perf.tau = parameters.tau;
                    truncated_performance(i) = performance_calculated(frec, fmin, options_perf);
                end
            end
            performance = max(truncated_performance);
            % performance = (truncated_performance(2) - truncated_performance(1))...
            %     * (value(6) - floor(value(6))) + truncated_performance(1);
        end
    end
end

if strcmpi(parameters.solvers_name(1), "rbds")
    if value(1) < 1 || (value(2) <= 0 && value(2) >= 1) ...
            || (value(3) < 0 || value(3) > value(4) || value(3) > value(5))...
            || (value(4) <= 0 || value(4) > value(5)) || value(5) <= 0 ...
            || value(6) < 0
        performance = NaN;
    else
        options_perf.natural_stop = false;
        if isinteger(value(6))
            parameters.solvers_options{1}.replacement_delay = value(6);
            [frec, fmin] = hp_calculated(parameters);
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
        else
            truncated_performance = NaN(2, 1);
            for i = 1 : 2
                parameters.solvers_options{1}.replacement_delay = floor(value(6) - 1) + i;
                [frec, fmin] = hp_calculated(parameters);
                if length(parameters.tau) > 1
                    num_tau = length(parameters.tau);
                    multi_performance = NaN(num_tau, 1);
                    for j = 1 : num_tau
                        options_perf.tau = parameters.tau(i);
                        multi_performance(j) = performance_calculated(frec, fmin, options_perf);
                    end
                    truncated_performance(i) = max(multi_performance);
                else
                    options_perf.tau = parameters.tau;
                    truncated_performance(i) = performance_calculated(frec, fmin, options_perf);
                end
            end
            performance = max(truncated_performance);
            % performance = (truncated_performance(2) - truncated_performance(1))...
            %     * (value(6) - floor(value(6))) + truncated_performance(1);
        end
    end
end

end


