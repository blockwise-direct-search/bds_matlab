function [x,fval,exitflag,output] = matlab_fminunc(FUN, x0, options)
% Find minimum of unconstrained multivariable function (using BFGS method)
global fhist_fminunc
global xhist_fminunc
fhist_fminunc = [];
xhist_fminunc = [];
[x,fval,exitflag,output] = fminunc(@(x)obj_fminunc(FUN,x),x0, options);
output.fhist = fhist_fminunc;
output.xhist = xhist_fminunc;