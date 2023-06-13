options.maxfun = 1e6;
options.StepTolerance = eps;

fullpath = mfilename('fullpath');
[path_examples,~] = fileparts(fullpath);
[path_bds, ~, ~] = fileparts(path_examples);
path_src = fullfile(path_bds, 'src');
addpath(path_src)

p = macup('akiva');
options.StepTolerance = 1e-10;

[x, fval, exitflag, output] = blockwise_direct_search(p.objective, p.x0, options);

fhist_length = length(output.fhist);
g_hist = NaN(1,fhist_length);
for eval_g = 1:fhist_length
    [~,gradient] = p.objective(output.xhist(:,eval_g));
    g_hist(eval_g) = norm(gradient);
end
gval = min(g_hist);
ratio = abs(gval)/options.StepTolerance

rmpath(path_src)


