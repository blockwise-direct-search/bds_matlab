function  [] = comparison(function_name)
fullpath = mfilename('fullpath');
[path_examples,~] = fileparts(fullpath);
[path_bds, ~, ~] = fileparts(path_examples);
path_src = fullfile(path_bds, 'src');
path_tests_private = fullfile(path_bds, 'tests', 'private');
path_competitors = fullfile(path_bds, 'tests', 'competitors');
addpath(path_src)
addpath(path_competitors)

old_dir = pwd;
cd(path_tests_private)
% Tell matlab where to find prima.
locate_prima();
locate_matcutest();
cd(old_dir)

p = macup(function_name);
options.maxfun = min(1e5, 1e3*length(p.x0));
options.StepTolerance = 0;
options.rhoend = options.StepTolerance;
fun = @(x) (1+1e-6*randn(1))*p.objective(x);

[x, fval, exitflag, output] = rbds(p.objective, p.x0, options)
%[x, fval, exitflag, output] = newuoa(p.objective, p.x0, options)

if isfield(output, 'xhist')
    fhist_length = length(output.fhist);
    g_hist = NaN(1,fhist_length);
    for eval_g = 1:fhist_length
        [~,gradient] = p.objective(output.xhist(:,eval_g));
        g_hist(eval_g) = norm(gradient);
    end
    gval = min(g_hist);
    ratio = abs(gval)/options.StepTolerance
end

cd(old_dir)
rmpath(path_src)
rmpath(path_competitors)
