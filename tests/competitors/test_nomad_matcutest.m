clear all
options_s2mpj.problem_type = 'u';
options_s2mpj.mindim = 1;
options_s2mpj.maxdim = 5;
problem_names = secup(options_s2mpj);
for i = 1:length(problem_names)
    p = macup(char(problem_names(i)));
    fprintf('Solving %s\n', p.name);
    nomad_wrapper(p.objective, p.x0, struct());
end