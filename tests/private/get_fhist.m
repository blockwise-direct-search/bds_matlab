function fhist_perfprof = get_fhist(p, maxfun_frec, j, r, solvers_options, test_options)
% GET_FHIST gets return value of j-th solver on the r-th randomized 
% experiment of problem p.

options = solvers_options{j};
name_solver = options.solver;
solver = str2func(name_solver);
% Initialze fhist for performance profile.
fhist_perfprof = NaN(maxfun_frec, 1);

% Scaling_matrix
% Gradient will be affected by scaling_matrix
% Find better way to deal with scaling_matrix
% if test_options.scaling_matrix
%    scaling_matrix = get_scaling_matrix(p,test_options);
% else
%    scaling_matrix = eye(length(p.x0));
% end

% Set maxfun before computing.
% Set MAXFUN to the maximum number of function evaluations if there exist
% some related parameters input, which are maxfun_factor and maxfun.
n = length(p.x0);
if isfield(test_options, "maxfun_factor") && isfield(test_options, "maxfun")
    options.maxfun = min(test_options.maxfun_factor*n, test_options.maxfun);
elseif isfield(test_options, "maxfun_factor")
    options.maxfun = test_options.maxfun_factor*n;
elseif isfield(test_options, "maxfun")
    options.maxfun = test_options.maxfun;
end

% Get the options for j-th solver.
[options] = get_options(name_solver, options);

% Turn off warning to save computation resource.
prima_list = ["cobyla", "uobyqa", "newuoa", "bobyqa", "lincoa"];
if ~isempty(find(prima_list == name_solver, 1))
    warnoff(name_solver);
end

% Try ... catch is to avoid stopping by the collapse of solvers. When some
% solver fails, we will use the iterates before it to record the fhist.
obj = ScalarFunction(p);
% try
%     solver(@(x)obj.fun(x,test_options.is_noisy,r,test_options), p.x0, options);
% catch
% end
solver(@(x)obj.fun(x,test_options.is_noisy,r,test_options), p.x0, options);
% Turn off warning is a very dangerous thing. So it must be set a loop to
% trun on after ending the computation.
if ~isempty(find(prima_list == name_solver, 1))
    warnoff(name_solver);
end

% Get length of fhist.
fhist_length = obj.nEval;
fhist_perfprof(1:fhist_length) = obj.valHist(1:fhist_length);

% Trim fhist for performance profile. If the length of fhist is less than maxfun,
% then the prolonged parts will be imparted the value of the last function evaluation 
% that we get eventually.
if  fhist_length < maxfun_frec
    if fhist_length > 0
        fhist_perfprof(fhist_length+1:maxfun_frec) = fhist_perfprof(fhist_length);
    end
else
    fhist_perfprof = fhist_perfprof(1:maxfun_frec);
end

end
