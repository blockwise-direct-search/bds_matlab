function p = s2mpj_wrapper(problem_info)

p.objective = @(x) problem_info.fun(x);
p.x0 = problem_info.x0;
p.name = problem_info.name;

end

