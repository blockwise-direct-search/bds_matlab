function x_best = ds_xin(fun, x0, options)

    if nargin < 3
        options = struct();
    end
    
    n = length(x0);
    D = [eye(n) -eye(n)];
    alpha = 1;
    expand = 2;
    shrink = 0.5;
    x_best = x0;
    f_best = fun(x_best);

    while alpha > eps
        for i = 1:2*n
            d = D(:, i);
            f_tmp = fun(x_best + alpha * d);
            if f_best - f_tmp > 1e-3 * alpha^2
                f_best = f_tmp;
                x_best = x_best + alpha * d;
                alpha = alpha * expand;
            else
                alpha = alpha * shrink;
                if alpha <= eps
                    break
                end
            end
        end
    end

end