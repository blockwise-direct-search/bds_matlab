function [fhist_perfprof, fval] = get_fhist(p, maxfun, j, r, solver_options, test_options)
% GET_FHIST gets return value of j-th solver on the r-th randomized 
% experiment of problem p.

name_solver = solver_options.solvers(j);
solver = str2func(name_solver);
% Initialze fhist for performance profile.
fhist_perfprof = NaN(maxfun, 1);

% Scaling_matrix
% Gradient will be affected by scaling_matrix
% Find better way to deal with scaling_matrix
% if test_options.scaling_matrix
%    scaling_matrix = get_scaling_matrix(p,test_options);
% else
%    scaling_matrix = eye(length(p.x0));
% end

% Set maxfun before computing.
if isfield(solver_options, "maxfun_dim")
   options.maxfun = min(solver_options.maxfun, solver_options.maxfun_dim*length(p.x0));
end

% Get the necessary options for j-th solver.
[options] = get_options(p, j, name_solver, solver_options, options);

% Turn off warning to save computation resource.
prima_list = ["cobyla", "uobyqa", "newuoa", "bobyqa", "lincoa"];
if ~isempty(find(prima_list == name_solver, 1))
    warnoff(name_solver);
end

% Try ... catch is to avoid stopping by the collapse of solvers. When some
% solver fails, we will use the iterates before it to record the fhist.
obj = ScalarFunction(p);
try
    solver(@(x)obj.fun(x,test_options.is_noisy,r,test_options), p.x0, options);
catch
end

% Turn off warning is a very dangerous thing. So it must be set a loop to
% trun on after ending the computation.
if ~isempty(find(prima_list == name_solver, 1))
    warnoff(name_solver);
end

% Fval should be the minimum among history of function values. Also, fval
% should always be the one without noise!
fval = min(obj.valHist);

% Get length of fhist.
fhist_length = obj.nEval;
fhist_perfprof(1:fhist_length) = obj.valHist(1:fhist_length);

% Trim fhist for performance profile. If the length of fhist is less than maxfun,
% then the prolonged parts will be imparted the value of the last function evaluation 
% that we get eventually.
if  fhist_length < maxfun
    fhist_perfprof(fhist_length+1:maxfun) = fhist_perfprof(fhist_length);
else
    fhist_perfprof = fhist_perfprof(1:maxfun);
end

end
