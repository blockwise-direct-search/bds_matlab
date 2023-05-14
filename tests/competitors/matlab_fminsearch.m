function [x, fval, exitflag, output] = matlab_fminsearch(FUN, x0, options)
% Find minimum of unconstrained multivariable function using
% derivative-free method (simplex method)
global fhist_fminsearch
global xhist_fminsearch
fhist_fminsearch = [];
xhist_fminsearch = [];
[x,fval,exitflag,output] = fminsearch(@(x)obj_fminsearch(FUN,x),x0, options);
output.fhist = fhist_fminsearch;
output.xhist = xhist_fminsearch;