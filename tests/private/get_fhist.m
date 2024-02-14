function [fhist, fhist_perfprof] = get_fhist(p, MaxFunctionEvaluations_frec, j, r, solvers_options, test_options)
% GET_FHIST gets return value of j-th solver on the r-th randomized
% experiment of problem p.

options = solvers_options{j};
name_solver = options.solver;
solver = str2func(name_solver);
% Initialize fhist for performance profile.
fhist_perfprof = NaN(MaxFunctionEvaluations_frec, 1);
% Set with_gradient for test_options. If the solver is fminunc_wrapper, then
% and the problem is noisy, then we will set with_gradient to be true.
% Make sure that the value of with_gradient is consistent with the value of
% with_gradient in the get_options.m. The only case for which we set
% with_gradient to be true is when the solver is fminunc_wrapper and the
% problem is noisy. 
test_options.with_gradient = strcmpi(name_solver, "fminunc_wrapper") && test_options.is_noisy;

% Try ... catch is to avoid stopping by the collapse of solvers. When some
% solver fails, we will use the iterates before it to record the fhist.
obj = ScalarFunction(p);

% try 
%     solver(@(x)obj.fun(x,test_options.is_noisy,r,test_options), p.x0, options);
% catch ME
%     warning(ME.identifier, '%s', ME.message);
%     warning('!!!Solver %s RAISE AN ERROR on problem %s with r = %d!!!', name_solver, p.name, r);
% end
solver(@(x)obj.fun(x,test_options.is_noisy,r,test_options), p.x0, options);

fhist = obj.valHist;
% Get length of fhist.
fhist_length = length(obj.valHist);
fhist_perfprof(1:fhist_length) = fhist;

% Trim fhist for performance profile. If the length of fhist is less than MaxFunctionEvaluations,
% then the prolonged parts will be imparted the value of the last function evaluation
% that we get eventually.
if  fhist_length < MaxFunctionEvaluations_frec
    if fhist_length > 0
        fhist_perfprof(fhist_length+1:MaxFunctionEvaluations_frec) = fhist_perfprof(fhist_length);
    end
else
    fhist_perfprof = fhist_perfprof(1:MaxFunctionEvaluations_frec);
end

end