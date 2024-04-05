function [f, f_real] = eval_fun(fun, x)
%EVAL_FUN Evaluate the function at the given point and it returns the real function value f and
%   the moderate function value f_real. f is used to record the history of the function values.
%   f_real is used for the optimization algorithm to make decisions, which deals with NaN, huge values, 
%   and evaluation failures.

try
    f_real = fun(x);    
catch
    warning('The function evaluation failed.');
    f_real = NaN; 
end

% Apply the moderate extreme barrier to handle NaN, huge values, and evaluation failures.
% See 4.5 of "PDFO: A Cross-Platform Package for Powell's Derivative-Free Optimization Solvers" 
% by Tom M. Ragonneau and Zaikun Zhang.
if isnan(f_real)
    f_real = inf;
end
f = min([f_real, 2^100, sqrt(realmax())]);

end

