function plot_fhist(problem_name, parameters)
% This file is to draw the function value history of the solvers.
% 

fullpath = mfilename("fullpath");
path_examples = fileparts(fullpath);
path_bds = fileparts(path_examples);
path_src = fullfile(path_bds, "src");
path_competitors = fullfile(path_bds, "tests", "competitors");
addpath(path_src)
addpath(path_competitors)

locate_matcutest();
p = macup(problem_name);

color_set = ["red", "blue", "green", "yellow"];
solvers_num = length(parameters.solvers_name);
fhist = cell(1, solvers_num);

if ~isfield(parameters, "solvers_options")
    parameters.solvers_options = {};
end

for i = 1:solvers_num
    parameters.solvers_options{i}.solver = parameters.solvers_name(i);
end

parameters = get_solvers(parameters);

for i = 1:solvers_num
    
    solver = str2func(parameters.solvers_options{i}.solver);
    if i <= length(parameters.solvers_options)
        options = parameters.solvers_options{i};
    end

    if strcmpi(parameters.solvers_name(i), "bds")
        [~, ~, ~, output] = solver(p.objective, p.x0, options);
    else
        obj = ScalarFunction(p);
        solver(@(x)obj.fun(x,false,0,struct()), p.x0);
        output.fhist = obj.valHist;
    end  
 
    fhist{i} = output.fhist;
end

for i = 1:solvers_num
    loglog(fhist{i}, color_set(i));
    hold on
end

rmpath(path_src)
rmpath(path_competitors)

end

