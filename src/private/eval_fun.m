function [f] = eval_fun(fun, x)
%EVAL_FUN evaluates function FUN at point X. If FUN is not well defined at X, return NaN. 
%

try
    f = fun(x);
catch
    % TODO: Is this the best value?
    f = NaN;  
end

end

