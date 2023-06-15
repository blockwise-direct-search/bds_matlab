function  [] = comparison(function_name)
fullpath = mfilename('fullpath');
[path_examples,~] = fileparts(fullpath);
[path_bds, ~, ~] = fileparts(path_examples);
path_src = fullfile(path_bds, 'src');
path_tests_private = fullfile(path_bds, 'tests', 'private');
path_competitors = fullfile(path_bds, 'tests', 'competitors');
addpath(path_src)
addpath(path_competitors)

path_temporary = pwd;
cd(path_tests_private)
% Tell matlab where to find prima.
locate_prima();
cd(pwd)

p = macup(function_name);
options.maxfun = min(1e5, 1e3*length(p.x0));
options.StepTolerance = 1e-10;
options.rhoend = options.StepTolerance;
fun = @(x) (1+1e-6*randn(1))*p.objective(x);

[x, fval, exitflag, output] = matlab_fminsearch(fun, p.x0, options)
[x, fval, exitflag, output] = newuoa(fun, p.x0, options)

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

rmpath(path_src)
rmpath(path_competitors)
