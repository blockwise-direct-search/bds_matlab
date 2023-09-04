function [fval,x] = obj_patternsearch(f,x)
global fhist_patternsearch
global xhist_patternsearch
xhist_patternsearch(:,end+1) = x; %Store the current trial point
fval = f(x);
fhist_patternsearch(end+1) = fval; %Store the current function evaluation
end

