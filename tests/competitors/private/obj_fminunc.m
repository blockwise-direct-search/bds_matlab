function [fval,x] = obj_fminunc(f,x)
global fhist_fminunc
global xhist_fminunc
xhist_fminunc(:,end+1) = x; %Store the current trial point
fval = f(x);
fhist_fminunc(end+1) = fval; %Store the current function evaluation
end

