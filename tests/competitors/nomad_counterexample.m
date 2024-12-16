function x = nomad_counterexample()

%wrapped_obj = @(varargin)objective(varargin{:});

x0 = [2,1];
fun = @(varargin) objective(varargin{:})

[x, ~, ~, ~, ~] = nomadOpt(fun, x0, -inf(2,1), inf(2,1), struct('MAX_BB_EVAL', '500', 'max_eval', '500'));

end



function f = objective(varargin)
    x = varargin{1};
    f = norm(x);
end


