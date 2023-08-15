fullpath = mfilename('fullpath');
[path_examples,~] = fileparts(fullpath);
[path_bds, ~, ~] = fileparts(path_examples);
path_src = fullfile(path_bds, 'src');
path_competitors = fullfile(path_bds, 'tests', 'competitors');
addpath(path_src)
addpath(path_competitors)

p = macup('akiva');
% p = macup('PALMER5C');
% p = macup('HEART6LS');
% p = macup('LANCZOS1LS');
options.StepTolerance = 1e-10;
options.Algorithm = "cbds";

[x, fval, exitflag, output] = bds(p.objective, p.x0, options);

fhist_length = length(output.fhist);
g_hist = NaN(1,fhist_length);
for eval_g = 1:fhist_length
    [~,gradient] = p.objective(output.xhist(:,eval_g));
    g_hist(eval_g) = norm(gradient);
end
gval = min(g_hist);
ratio = abs(gval)/options.StepTolerance

rmpath(path_src)
rmpath(path_competitors)
