function matlab_fminsearch(FUN, x0, options)
% Find minimum of unconstrained multivariable function using
% derivative-free method (simplex method)
%

fminsearch(@(x)obj_fminsearch(FUN,x),x0, options);
