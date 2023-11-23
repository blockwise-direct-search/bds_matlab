function [f] = eval_fun(fun, x)
%EVAL_FUN evaluates function FUN at point X. If FUN is not well defined at X, return NaN. 
%

try
    f = fun(x);    
catch
    warning('The function evaluation failed.');
    f = NaN; 
end

% Apply the moderate extreme barrier to handle NaN, huge values, and evaluation failures.
% See 4.5 of "PDFO: A Cross-Platform Package for Powell's Derivative-Free Optimization Solvers" 
% by Tom M. Ragonneau and Zaikun Zhang.
if isnan(f)
    f = inf;
end
f = min([f, 2^100, sqrt(realmax())]);

end

