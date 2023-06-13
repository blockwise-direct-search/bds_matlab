options.maxfun = 1e4;
options.StepTolerance = eps;

fullpath = mfilename('fullpath');
[path_examples,~] = fileparts(fullpath);
[path_bds, ~, ~] = fileparts(path_examples);
path_src = fullfile(path_bds, 'src');
add(path_src)

[x, fval, exitflag, output] = blockwise_direct_search(@rosenb, [0; 0; 0], options)

rmpath(path_src)

function f = rosenb(x)

f = 0;
for i = 1:(length(x)-1)
	f = f+100*(x(i+1)-x(i)^2)^2+(x(i)-1)^2;
end

end

