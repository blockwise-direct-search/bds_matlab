clear all
fullpath = mfilename('fullpath');
[path_examples,~] = fileparts(fullpath);
[path_bds, ~, ~] = fileparts(path_examples);
path_src = fullfile(path_bds, 'src');
path_competitors = fullfile(path_bds, 'tests', 'competitors');
addpath(path_src)
addpath(path_competitors)

p = macup('akiva');
random_seed = abs(ceil(1e5*sin(sum(p.x0)))) + ...
    abs(ceil(1e4 * sin(1e3*norm(p.x0)))) + 5000 * norm(p.x0);
noise_level = 1e-3;
with_failure = true;

options.StepTolerance = 1e-10;
options.Algorithm = "cbds";

tough_problem = tough(p, random_seed, noise_level, with_failure);

[x, fval, exitflag, output] = bds(tough_problem.objective, tough_problem.x0, options);

rmpath(path_src)
rmpath(path_competitors)
