function get_fhist_one_problem(p, solvers_options, test_options)

if nargin < 2
    solvers_options = {};
    test_options = struct();
end

hfig = figure("visible", true);
MaxFunctionEvaluations_frec = 500*length(p.x0);

% Set up the test options.
if ~isfield(test_options, "is_noisy")
    test_options.is_noisy = true;
end
if ~isfield(test_options, "noise_level")
    test_options.noise_level = 1e-4;
end
test_options.is_abs_noise = false;
test_options.noise_type = "gaussian";
if ~isfield(test_options, "num_random")
    test_options.num_random = 1;
end

if isempty(solvers_options)
    solvers_options{1}.solver = "bds";
    solvers_options{1}.Algorithm = "ds";

    solvers_options{2}.solver = "prima_wrapper";
    solvers_options{2}.Algorithm = "newuoa";
    solvers_options{2}.MaxFunctionEvaluations = MaxFunctionEvaluations_frec;
end

%solvers_options{2}.fminunc_type = "bfgs";
%solvers_options{2}.with_gradient = test_options.is_noisy;

% Set up the function handle. If p is a Problem object, then we neet to convert it to a struct.
if isa(p, 'Problem')
    problem_orig = str2func(char(p.name));
    problem_info = problem_orig('setup');
    p = s2mpj_wrapper(problem_info, p.name);
end

[fhist_bds, ~] = get_fhist(p, MaxFunctionEvaluations_frec, 1, 1, solvers_options, test_options);
[fhist_newuoa, ~] = get_fhist(p, MaxFunctionEvaluations_frec, 2, 1, solvers_options, test_options);

% Get fval.
fval = min(min(fhist_newuoa), min(fhist_bds));

% Deal fhist.
% cummin is used to get the best function value found so far.
fhist_bds_plot = cummin(fhist_bds);
fhist_newuoa_plot = cummin(fhist_newuoa);

fhist_bds_plot = abs(fhist_bds_plot - fval + eps)/max(abs(fhist_bds(1)), eps);
fhist_newuoa_plot = abs(fhist_newuoa_plot - fval + eps)/max(abs(fhist_newuoa(1)), eps);


% Plot the figure.
opengl('OpenGLSoftware');
set(gcf, 'Color', 'w'); % Set background color to white.
semilogy(fhist_bds_plot, 'blue');
hold on;
semilogy(fhist_newuoa_plot, 'red');
hold off;

% Add legend, title, and labels.
solvers_legend = {'ds', 'newuoa'};
legend(solvers_legend{1}, solvers_legend{2});
if test_options.is_noisy
    title_name = strcat(p.name, "-", "noisy", num2str(log10(test_options.noise_level)));
else
    title_name = p.name;
end
title(title_name);
xlabel('Function evaluations');
ylabel('Function values');

% Save the figure.
currentFilePath = mfilename('fullpath');
path_testdata = fullfile(fileparts(fileparts(currentFilePath)), 'testdata');
keyboard
if test_options.is_noisy
    eps_name = strcat(p.name, "_", solvers_legend{1}, "_", solvers_legend{2}, '_noisy', "_", num2str(log10(test_options.noise_level)), '.eps');
else
    eps_name = strcat(p.name, "_", solvers_legend{1}, "_", solvers_legend{2}, '_noiseless', '.eps');
end
saved_path = fullfile(path_testdata, eps_name);
exportgraphics(gcf, saved_path, 'ContentType', 'vector');

end

