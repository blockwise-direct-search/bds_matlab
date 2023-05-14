options.maxfun = 1e4;
options.tol = eps;

[x, fval, exitflag, output] = blockwise_direct_search(@rosenb, [0; 0; 0], options)

function f = rosenb(x)

f = 0;
for i = 1:(length(x)-1)
	f = f+100*(x(i+1)-x(i)^2)^2+(x(i)-1)^2;
end

end

