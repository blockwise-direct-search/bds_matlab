% This script is for tough test.

% Record the current path.
oldpath = path(); 
% Restore the "right out of the box" path of MATLAB. 
restoredefaultpath;  
% Record the current directory.
old_dir = pwd();

fullpath = mfilename('fullpath');
[path_examples,~] = fileparts(fullpath);
[path_bds, ~, ~] = fileparts(path_examples);
path_src = fullfile(path_bds, 'src');
path_competitors = fullfile(path_bds, 'tests', 'competitors');
addpath(path_src);
addpath(path_competitors);

% Tell MATLAB where to find MatCUTEst.
locate_matcutest();

% Get list of problems
if parameters.problems_dim == "small"
    s.type = "u"; % Unconstrained: 'u'
    s.mindim = 1; % Minimum of dimension
    s.maxdim = 5; % Maximum of dimension
elseif parameters.problems_dim == "big"
    s.type = "u"; % Unconstrained: 'u'
    s.mindim = 6; % Minimum of dimension
    s.maxdim = 100; % Maximum of dimension
end

s.blacklist = [];
% Problems that crash.
s.blacklist = [s.blacklist, {}];
s.blacklist = [{'LRCOVTYPE'},{'LRIJCNN1'},{'PARKCH'},{'STRATEC'}];

problem_names = secup(s);
num_problems = length(problem_names);
fprintf("We will load %d problems\n\n", num_problems);

options.Algorithm = parameters.Algorithm;
noise_level = 1e-3;

for i = 1:num_problems
    p = macup(problem_names(1, i));
    fprintf("%s\n", p.name);
    random_seed = abs(ceil(1e5*sin(sum(p.x0)))) + ...
        abs(ceil(1e4 * sin(1e3*norm(p.x0)))) + 5000 * norm(p.x0);
    with_failure = true;
    tough_problem = tough(p, random_seed, noise_level, with_failure);
    [x, fval, exitflag, output] = bds(tough_problem.objective, tough_problem.x0, options);
end

% Restore the path to oldpath.
setpath(oldpath);  
% Go back to the original directory.
cd(old_dir);

