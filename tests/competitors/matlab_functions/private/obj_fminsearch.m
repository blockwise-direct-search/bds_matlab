function [fval,x] = obj_fminsearch(f,x)
global fhist_fminsearch
global xhist_fminsearch
xhist_fminsearch(:,end+1) = x; %Store the current trial point
fval = f(x);
fhist_fminsearch(end+1) = fval; %Store the current function evaluation
end

