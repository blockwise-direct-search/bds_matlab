function [xopt, fopt, exitflag, output] = test_problem(Algorithm)
%Test bds on some specific problems.
%

if nargin < 1
    options = struct();
end

options.Algorithm = Algorithm;
options.maxfun = 1e5;
options.StepTolerance = eps;
options.output_xhist = true;
options.output_alpha_hist = true;
options.output_block_hist = true;
[xopt, fopt, exitflag, output] = bds(@(x)hmlb(x), [0;0], options)
end

function f = goldp(x)
%GOLDP evaluates the Goldstein-Price function
%
%   See
%   [1] Dixon, L. C. W., & Szego, G. P. (1978). The global optimization problem: an introduction. Towards global optimization, 2, 1-15.

f1a = (x(1) + x(2) + 1)^2;
f1b = 19 - 14*x(1) + 3*x(1)^2 - 14*x(2) + 6*x(1)*x(2) + 3*x(2)^2;
f1 = 1 + f1a*f1b;

f2a = (2*x(1) - 3*x(2))^2;
f2b = 18 - 32*x(1) + 12*x(1)^2 + 48*x(2) - 36*x(1)*x(2) + 27*x(2)^2;
f2 = 30 + f2a*f2b;

f = f1*f2;

return
end

function f = chebquad(x)
%CHEBQUAD evaluates the Chebyquad function.
%
%   See
%   [1] Fletcher (1965), "Function minimization without evaluating derivatives --- a review"

n = length(x);
y(1,1:n) = 1;
y(2, 1:n) = 2*x(1:n) - 1;
for i = 2:n
    y(i+1, 1:n) = 2*y(2, 1:n).*y(i, 1:n) - y(i-1, 1:n);
end
f = 0;
for i = 1 : n+1
    tmp = mean(y(i, 1:n));
    if (mod(i, 2) == 1)
        tmp=tmp+1/double(i*i-2*i);
    end
    f = f + tmp*tmp;
end

return
end

function [f, g] = hmlb(x)
%HMLB evaluates the Himmelblau"s function and its gradient
%
%   See
%   [1]  Himmelblau (1972),  "Applied Nonlinear Programming"

f = (x(1)^2+x(2)-11)^2 + (x(1)+x(2)^2-7)^2;
g = 2*[-7 + x(1) + x(2)^2 + 2*x(1)*(-11 + x(1)^2 + x(2)); -11 + x(1)^2 + x(2) + 2*x(2)*(-7 + x(1) + x(2)^2)];

return
end

function f = mcc(x)
%MCC evaluates the McCormick function

f1 = sin(x(1) + x(2));
f2 = (x(1) - x(2))^2;
f3 = -1.5*x(1);
f4 = 2.5*x(2);

f = f1 + f2 + f3 + f4 + 1;

return
end