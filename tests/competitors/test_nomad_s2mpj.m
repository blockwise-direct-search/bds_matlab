clear all
options_s2mpj.problem_type = 'u';
options_s2mpj.mindim = 1;
options_s2mpj.maxdim = 5;
problem_names = s_select(options_s2mpj);
for i = 1:length(problem_names)
    q = s_load(char(problem_names(i)));
    fprintf('Solving %s\n', q.name);
    p.objective = @(x) q.fun(x);
    p.x0 = q.x0;
    keyboard
    nomad_wrapper(p.objective, p.x0, struct());
end