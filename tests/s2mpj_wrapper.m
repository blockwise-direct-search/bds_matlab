function problem_struct = s2mpj_wrapper(q, problem_name)

problem_name_fun = str2func(char(problem_name));
problem_struct = struct();
problem_struct.name = char(problem_name);
problem_struct.x0 = q.x0;
problem_struct.objective = @(x) problem_name_fun('fx', x);

end

