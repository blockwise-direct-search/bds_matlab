function get_fhist_one_problem(p, solvers_options, test_options, saved_dir)

if nargin < 2
    solvers_options = {};
    test_options = struct();
    saved_path = '';
end


hfig = figure("visible", true);
MaxFunctionEvaluations_frec = 500*length(p.x0);

% Set up the test options.
if ~isfield(test_options, "is_noisy")
    test_options.is_noisy = true;
end
if ~isfield(test_options, "noise_level")
    test_options.noise_level = 1e-1;
end
test_options.is_abs_noise = false;
test_options.noise_type = "gaussian";
if ~isfield(test_options, "num_random")
    test_options.num_random = 2;
end

if isempty(solvers_options)
    solvers_options{1}.solver = "bds";
    solvers_options{1}.Algorithm = "ds";

    solvers_options{2}.solver = "prima_wrapper";
    solvers_options{2}.Algorithm = "newuoa";
    solvers_options{2}.MaxFunctionEvaluations = MaxFunctionEvaluations_frec;
    solvers_options{2}.iprint = 2;
end

%solvers_options{2}.fminunc_type = "bfgs";
%solvers_options{2}.with_gradient = test_options.is_noisy;

% Set up the function handle. If p is a Problem object, then we neet to convert it to a struct.
if isa(p, 'Problem')
    problem_orig = str2func(char(p.name));
    problem_info = problem_orig('setup');
    p = s2mpj_wrapper(problem_info, p.name);
end

[fhist_bds, ~] = get_fhist(p, MaxFunctionEvaluations_frec, 1, test_options.num_random, solvers_options, test_options);
[fhist_newuoa, ~] = get_fhist(p, MaxFunctionEvaluations_frec, 2, test_options.num_random, solvers_options, test_options);

% Get fval.
fval = min(min(fhist_newuoa), min(fhist_bds));

% Deal fhist.
% Shift the function values to make them positive.
shiftValue = abs(fval) + test_options.noise_level;
fhist_bds_plot = fhist_bds + shiftValue;
fhist_newuoa_plot = fhist_newuoa + shiftValue;

% Plot the figure.
opengl('OpenGLSoftware');
set(gcf, 'Color', 'w'); % Set background color to white.
semilogy(fhist_bds_plot, 'blue');
hold on;
semilogy(fhist_newuoa_plot, 'red');
hold off;

% Restore y-axis ticks to original values
y_ticks_shifted = get(gca, 'YTick'); % Get current y-axis ticks
set(gca, 'YTick', y_ticks_shifted - shiftValue); % Shift back to original values

% Set the y-axis ticks to be in the order of 10^k.
y_lim = ylim; % Get the y-axis limits
y_upper = ceil(log10(y_lim(2))); % Obtain the upper limit of the y-axis
y_lower = floor(log10(y_lim(1)));  % Obtain the lower limit of the y-axis

% Establish the new y-axis limits.
new_ylim = [10^(y_lower), 10^(y_upper)];
% Set the y-axis limits to be the new limits.
ylim(new_ylim);

% Set the distance between the ticks to be 2.
if y_upper - y_lower <= 1
    selected_ticks = 10.^(y_lower:1:y_upper);
else
    selected_ticks = 10.^(y_lower:2:y_upper);
end

% Set the y-axis ticks to be the selected ticks.
set(gca, 'YTick', selected_ticks);

% Adjust font size as needed
set(gca, 'FontSize', 10);
grid off;

% Set x-axis ticks
scale = length(p.x0) + 1;
x_max = max(length(fhist_bds_plot), length(fhist_newuoa_plot)); 
x_ticks = 0:floor(x_max/5):x_max;

set(gca, 'XTick', x_ticks); 
set(gca, 'XTickLabel', floor(x_ticks / scale));
xlim([0 x_max]);


% Add legend, title, and labels.
solvers_legend = {'ds', 'newuoa'};
legend(solvers_legend{1}, solvers_legend{2});
if test_options.is_noisy
    title_name = strcat(p.name, "-", num2str(length(p.x0)), "-", "noisy", num2str(log10(test_options.noise_level)));
else
    title_name = strcat(p.name, "-", num2str(length(p.x0)));
end
title(title_name);
xlabel('Function evaluations');
ylabel('Simplex gradient');

% Save the figure.
currentFilePath = mfilename('fullpath');
path_testdata = fullfile(fileparts(fileparts(currentFilePath)), 'testdata');
if test_options.is_noisy
    eps_name = strcat(p.name, "_", solvers_legend{1}, "_", solvers_legend{2}, '_noisy', "_", num2str(log10(test_options.noise_level)), '.eps');
else
    eps_name = strcat(p.name, "_", solvers_legend{1}, "_", solvers_legend{2}, '_plain', '.eps');
end

if ~isempty(saved_dir)
    saved_path = fullfile(saved_dir, eps_name);
else
    saved_path = fullfile(path_testdata, eps_name);
end

% Make sure the directory exists.
fig = gcf; % Get the current figure
% Export the figure to the saved path.
print(fig, saved_path, '-depsc'); % Save the figure as an eps file.

% % Hide the toolbar.
% set(gca, 'Toolbar', 'none');
% exportgraphics(gcf, saved_path, 'ContentType', 'vector');

end

% cummin is used to get the best function value found so far.
% fhist_bds_plot = cummin(fhist_bds);
% fhist_newuoa_plot = cummin(fhist_newuoa);

% fhist_bds_plot = abs(fhist_bds_plot - fval + eps)/max(abs(fhist_bds(1)), eps);
% fhist_newuoa_plot = abs(fhist_newuoa_plot - fval + eps)/max(abs(fhist_newuoa(1)), eps);