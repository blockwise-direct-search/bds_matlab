function [f, f_real] = eval_fun(fun, x)
%EVAL_FUN evaluates function FUN at point X, returning f and f_real.
%   f_real is the real function value, while f is a moderated version of f_real.
%   The moderation is to handle NaN, huge values, and evaluation failures. The 
%   algorithm will operate on f, while f_real is used for recording the history.

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
f = min([f_real, 10^30, sqrt(realmax())]);

end

