function [ratio, gval] = test_gradient_stepsize(problem_name, options)
% This function tests the ratio between the gradient norm and the StepTolerance on CUTEst problems.
% 

format long

if nargin < 2
    options.solver_name = "bds";
end

fullpath = mfilename('fullpath');
path_examples = fileparts(fullpath);
path_bds = fileparts(path_examples);
path_src = fullfile(path_bds, 'src');
path_locate = fullfile(path_bds, 'tests', 'private');
path_competitors = fullfile(path_bds, 'tests', 'competitors');
addpath(path_src)
addpath(path_competitors)

cd(path_locate)
locate_matcutest();
p = macup(problem_name);
if strcmpi(options.solver_name, "newuoa")
    locate_prima();
end

solver = str2func(options.solver_name);

tic;
if strcmpi(options.solver_name, "bds") || strcmpi(options.solver_name, "bfo_optimize") ||...
        strcmpi(options.solvers_name, "newuoa")
    [~, ~, ~, output] = solver(p.objective, p.x0, options);
else
    obj = ScalarFunction(p);
    solver(@(x)obj.fun(x,false,0,struct()), p.x0);
    output.fhist = obj.valHist;
end
toc;

fhist_length = length(output.fhist);
g_hist = NaN(1,fhist_length);
if isfield(output, "xhist")
    for eval_g = 1:fhist_length
        [~,gradient] = p.objective(output.xhist(:,eval_g));
        g_hist(eval_g) = norm(gradient);
    end
    gval = min(g_hist);
    ratio = abs(gval)/options.StepTolerance;
end

rmpath(path_src)
rmpath(path_competitors)

end


