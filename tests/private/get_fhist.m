function [fhist_perfprof, fval] = get_fhist(p, maxfun, j, r, options_solvers, options_test)
% Get fhist and related information of j-th solver on p problem.

% Gradient will be affected by scaling_matrix
name_solver = options_solvers.solvers(j);
solver = str2func(name_solver);
fhist_perfprof = NaN(maxfun,1);

% Scaling_matrix
% Find better way to deal with scaling_matrix
% if options_test.scaling_matrix
%    scaling_matrix = get_scaling_matrix(p,options_test);
% else
%    scaling_matrix = eye(length(p.x0));
% end

% Maxfun_dim will be an input.
if isfield(options_solvers, "maxfun_dim")
   options.maxfun = min(options_solvers.maxfun, options_solvers.maxfun_dim*length(p.x0));
end

[options] = get_options(p, j, name_solver, options_solvers, options);

% Turn off warning to save computation resource.
prima_list = ["cobyla", "uobyqa", "newuoa", "bobyqa", "lincoa"];
if ~isempty(find(prima_list == name_solver, 1))
    warnoff(name_solver);
end

% Experimence with noise (if num_random == 1, then the experiment has no noise)
obj = ScalarFunction(p);
solver(@(x)obj.fun(x,options_test.is_noisy,r,options_test),p.x0, options);

% Turn off warning is a very dangerous thing. So it must be set a loop to
% trun on after ending the computation.
if ~isempty(find(prima_list == name_solver, 1))
    warnoff(name_solver);
end

% In case meeting simple decrease but not sufficient decrease. Also, fval
% should always be the one without noise!
fval = min(obj.valHist);

% length of fhist and ghist
fhist_length = obj.nEval; 
% if isfield(output, "xhist")
%     g_hist = NaN(1,fhist_length);
%     for eval_g = 1:fhist_length
%         [~,gradient] = p.objective(output.xhist(:,eval_g));
%         g_hist(eval_g) = norm(gradient);
%     end
%     gval = min(g_hist);
%     gval_relative = gval/g_hist(1);
% end
fhist_perfprof(1:fhist_length) = obj.valHist(1:fhist_length);
if  fhist_length < maxfun
    fhist_perfprof(fhist_length+1:maxfun) = fhist_perfprof(fhist_length);
else
    fhist_perfprof = output_perfprof.fhist(1:maxfun);
end
end

% function FUN = objective(x, p, ~, ~, ~)
% FUN = p.objective(x);
% end 


% function FUN = objective_noise(x, p, k, scaling_matrix, options_test)
% if ~options_test.scaling_matrix
%     FUN = p.objective(x);
% else
%     FUN = p.objective(scaling_matrix*x);
% end
% % What is is_noisy, noise_type, noise_level, seed?
% if options_test.is_noisy
%     seed = length(x) + abs(ceil(1e4 * sin(k))) + 5000 * k;
%     rng(seed)
%     if ~strcmpi(options_test.noise_type, 'uniform')
%        noise = rand(1); 
%     end
%     if ~strcmpi(options_test.noise_type, 'gaussian')
%        noise = randn(1); 
%     end
%     rng(seed)
%     if options_test.noise_abs == "relative"
%         FUN = FUN*(1.0+options_test.noise_level*noise);
%     else
%         FUN = FUN + options_test.noise_level*noise;
%     end
% end
% return
% end 

