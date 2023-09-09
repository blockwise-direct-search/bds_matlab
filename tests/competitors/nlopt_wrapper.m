function [xval, fval, exitflag, output] = nlopt_wrapper(fun, x0, options)

global nlopt_fhist
global nlopt_nf

% Dimension
n = numel(x0);

%opt.algorithm = NLOPT_LN_COBYLA;
opt.algorithm = options.Algorithm;
opt.min_objective = @(x)objective(fun,x);

if ~isfield(options, 'maxfun')
    options.maxfun = 1000*n;
end

nlopt_nf = 0;
nlopt_fhist = NaN(1, options.maxfun);

[xval, ~, ~] = nlopt_optimize(opt, x0');
xval = xval';

fval = fun(xval);

exitflag = -1;

output.nf = nlopt_nf;
output.fhist = nlopt_fhist(1:nlopt_nf);

end

function [f] = objective(fun, x)

global nlopt_fhist
global nlopt_nf
f = fun(x);
nlopt_nf = nlopt_nf+1;
nlopt_fhist(nlopt_nf) = f;

end
