function matlab_patternsearch(FUN, x0, options)
% Find minimum of unconstrained multivariable function (using patternsearch
% method).
%

patternsearch(@(x)obj_patternsearch(FUN,x),x0,[],[],[],[],[],[],[],options);