function matlab_fminunc(FUN, x0, options)
% Find minimum of unconstrained multivariable function (using BFGS method).
%

fminunc(@(x)obj_fminunc(FUN,x),x0, options);
