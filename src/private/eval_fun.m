function [f] = eval_fun(fun, x)
%EVAL_FUN evaluates function FUN at point X. If FUN is not well defined at X, return NaN. 
%

try
    f = fun(x);    
catch
    % TODO: if there is a failure of the function, give a warning.
    f = NaN; 
end

% Apply the moderate extreme barrier to handle NaN, huge values, and evaluation failures.
% See 4.5 of "PDFO: A Cross-Platform Package for Powell's Derivative-Free Optimization Solvers" 
% by Tom M. Ragonneau and Zaikun Zhang.

q = quantizer('double');
maxexponent = exponentmax(q);
radix = 2;

if isnan(f)
    f = inf;
end

f = min(f, radix^min(100, maxexponent / 2));

end

