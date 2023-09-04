function fval = obj(f,x)
global fhi
fval = f(x);
fhi(end+1) = fval; %Store the current function evaluation
end

