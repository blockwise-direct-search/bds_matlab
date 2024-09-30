function performance_return = hp_handle(value, parameters)
% perfprof_handle - Function handle for calculating the value of the
% performance profile

if strcmpi(parameters.tuning_solver, "newuoa")
    if strcmpi(parameters.solvers_name(1), "pbds") || strcmpi(parameters.solvers_name(1), "rbds")
        value = [value(1:2); exp(value(3:end-1)) + eps; value(end)];
    else
        value = [value(1:2); exp(value(3:end)) + eps];
    end
end

if isfield(parameters, "reduction_error") && parameters.reduction_error
    performance_list = ones(7, 1);
else
    performance_list = 1;
end

% We only do central differences for the hyperparameters that are far smaller
% than the others, which are reduction_factors.
coordinate_matrix = eye(5);
for j = 1:length(performance_list)
    if j == 2 || j == 3
        alpha = abs(value(3)/2);
        hyperparameters_value = value + (-2*j + 5)*alpha*coordinate_matrix(:, 3);
    elseif j == 4 || j == 5
        alpha = abs(value(4)/2);
        hyperparameters_value = value + (-2*j + 9)*alpha*coordinate_matrix(:, 4);
    elseif j == 6 || j == 7
        alpha = abs(value(5)/2);
        hyperparameters_value = value + (-2*j + 13)*alpha*coordinate_matrix(:, 5);
    else
        hyperparameters_value = value;
    end

    penalty = 1000;
    dist = hp_regularization(hyperparameters_value);
    x_projected = hp_projection(hyperparameters_value);

    % Notice that the value of the parameters may not be the same as the ones used
    % in the tuning_solver. For example, log(exp(eps)) ~= eps numerically.
    parameters.solvers_options{1}.expand = x_projected(1);
    parameters.solvers_options{1}.shrink = x_projected(2);
    parameters.solvers_options{1}.reduction_factor = x_projected(3:5);
    
    if strcmpi(parameters.solvers_name(1), "pbds")
        parameters.solvers_options{1}.permuting_period = x_projected(end);
    end

    if strcmpi(parameters.solvers_name(1), "rbds")
        parameters.solvers_options{1}.replacement_delay = x_projected(end);
    end
    
    for i = 1:length(parameters.solvers_name)
        parameters.solvers_options{i}.solver = parameters.solvers_name(i);
    end

    [~, frec, fmin] = profile(parameters);
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

    performance_list(j) = performance;

end

performance_return = mean(performance_list);

end


