function [f] = eval_fun(fun, x)
%EVAL_FUN Evaluate function fun at point x. If the function is not well defined at x, return NaN. 
try
    f = fun(x);
catch
    % Is this the best value?
    f = NaN;  
end

end

