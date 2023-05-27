function [x,fval,exitflag,output] = matlab_patternsearch(FUN, x0, options)
% Find minimum of unconstrained multivariable function (using patternsearch method)
global fhist_patternsearch
global xhist_patternsearch
fhist_patternsearch = [];
xhist_patternsearch = [];
[x,fval,exitflag,output] = patternsearch(@(x)obj_patternsearch(FUN,x),x0,[],[],[],[],[],[],[],options);
output.fhist = fhist_patternsearch;
output.xhist = xhist_patternsearch;