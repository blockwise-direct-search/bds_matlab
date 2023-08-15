function [f] = eval_fun(fun, x)
%EVAL_FUN Summary of this function goes here
%   Detailed explanation goes here
try
    f = fun(x);
catch
    f = NaN;  % Is this the best value?
end
% Do you need to do something else?
end

