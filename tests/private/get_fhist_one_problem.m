function get_fhist_one_problem(p)

% The problem list is {'EG2'}, {'EXTROSNB'}, {'FLETCHCR'}, {'GENROSE'}, {'INTEQNELS'}, {'KSSLS'}, 
% {'OSCIPATH'}, {'PENALTY1'}, from the CUTEst test set (unconstrained problems, 500\leq n \leq 1000).
hfig=figure("visible", true);
MaxFunctionEvaluations_frec = 500*length(p.x0);

test_options.is_noisy = false;
test_options.noise_level = 1e-3;
test_options.is_abs_noise = false;
test_options.noise_type = "gaussian";
test_options.num_random = 1;

solvers_options{1}.solver = "bds";
solvers_options{1}.Algorithm = "cbds";

solvers_options{2}.solver = "fminunc_wrapper";
solvers_options{2}.fminunc_type = "bfgs";
solvers_options{2}.with_gradient = test_options.is_noisy;


[fhist_bds, ~] = get_fhist(p, MaxFunctionEvaluations_frec, 1, 1, solvers_options, test_options);
[fhist_bfgs, ~] = get_fhist(p, MaxFunctionEvaluations_frec, 2, 1, solvers_options, test_options);

% Get fval.
fval = min(min(fhist_bfgs), min(fhist_bds));

% Deal fhist.
fhist_bds_plot = cummin(fhist_bds);
fhist_bfgs_plot = cummin(fhist_bfgs);

fhist_bds_plot = abs(fhist_bds_plot - fval + eps)/max(abs(fhist_bds(1)), eps);
fhist_bfgs_plot = abs(fhist_bfgs_plot - fval + eps)/max(abs(fhist_bfgs(1)), eps);

semilogy(fhist_bds_plot, 'red');
hold on
semilogy(fhist_bfgs_plot, 'blue');
keyboard

end

