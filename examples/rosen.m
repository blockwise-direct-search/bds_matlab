options.maxfun = 1e4;
options.tol = eps;
addpath('/home/lhtian97/bds_new_framework/src');
[x, fval, exitflag, output] = blockwise_direct_search(@rosenb, [0; 0; 0], options)
rmpath('/home/lhtian97/bds_new_framework/src');

function f = rosenb(x)

f = 0;
for i = 1:(length(x)-1)
	f = f+100*(x(i+1)-x(i)^2)^2+(x(i)-1)^2;
end

end

