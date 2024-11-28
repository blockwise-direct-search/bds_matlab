function plot_parameters(parameters, solver, competitor, options)
% parameters: a structure with two fields; the field names are the names of the parameters; for each
% field, the value is a vector representing the values of the corresponding parameter.
% solver: a string representing the solver whose performance is to be evaluated.
% competitor: a string representing the competitor solver.
% options: a structure representing the options to be passed to the performance function.

% Get parameter names
param_names = fieldnames(parameters);
assert(length(param_names) == 2, 'There should be two parameters.');
param1_name = param_names{1};
param2_name = param_names{2};

% Create a grid of parameter values
[p1, p2] = meshgrid(parameters.(param1_name), parameters.(param2_name));

% Initialize performance matrix
perfs = NaN(size(p1));

% Get performance for each parameter combination
parfor ip = 1:numel(p1)
    % Set solver options
    solver_options = struct();
    solver_options.(param1_name) = p1(ip);
    solver_options.(param2_name) = p2(ip);

    % Pass solver_options to the performance function via local_options. The performance function
    % should then pass solver_options to the solver.
    local_options = options;
    local_options.solver_options = solver_options;

    % Compute performance
    fprintf('Evaluating performance for %s = %f, %s = %f\n', param1_name, p1(ip), param2_name, p2(ip));
    perfs(ip) = eval_performance(solver, competitor, local_options);
end


% Plot the save the results

% We save the results in the `data_path` folder. 
current_path = fileparts(mfilename("fullpath"));
% Create the folder if it does not exist.
data_path = fullfile(current_path, "tuning_data");
if ~exist(data_path, 'dir')
    mkdir(data_path);
end
% Creat a subfolder stamped with the current time for the current test. 
time_str = char(datetime('now', 'Format', 'yy_MM_dd_HH_mm'));
data_path = fullfile(data_path, time_str);
mkdir(data_path);

% Save performance data 
save(fullfile(data_path, 'performance_data.mat'), 'p1', 'p2', 'perfs');

% Plot
FigHandle=figure('Name', ['(', param1_name, ', ', param2_name, ')', ' v.s. performance']);
title(gca,'performance');
xlabel(param1_name);
ylabel(param2_name);
hold on;

colormap(jet);
if isfield(options, 'log_color') && options.log_color
    % Use log scale of perfs for a better usage of the color spectrum.
    max_perf = max(perfs(:));
    min_perf = min(perfs(:));
    C = min_perf + (max_perf - min_perf) .* log(perfs - min_perf + 1) ./ log(max_perf - min_perf + 1);
    surf(p1, p2, perfs, C, 'FaceColor','interp');
else
    surf(p1, p2, perfs, 'FaceColor','interp');
end
colorbar; 

view(3) % 3D view
% Save eps of 3d plot 
saveas(FigHandle, fullfile(data_path, [param1_name, '_', param2_name, '_vs_performance_3d.eps']), 'epsc');
% Try converting the eps to pdf.
epsPath = fullfile(data_path, [param1_name, '_', param2_name, '_vs_performance_3d.eps']);
system(('epstopdf '+epsPath));

% Save eps of 2d plot 
view(2); % Top-down view
saveas(FigHandle, fullfile(data_path, [param1_name, '_', param2_name, '_vs_performance_2d.eps']), 'epsc');
% Try converting the eps to pdf.
epsPath = fullfile(data_path, [param1_name, '_', param2_name, '_vs_performance_2d.eps']);
% One way to convert eps to pdf.
system(('epstopdf '+epsPath));

fprintf('Performance data and plots saved in \n %s\n', data_path);

end
