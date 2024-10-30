function [xopt, fopt, xhist, fhist, alpha_hist] = pds_xin(fun, x0, options)

    if nargin < 3
        options = struct();
    end

    if isfield(options, 'maxfun')
        maxfun = options.maxfun;
    else
        maxfun = 1000 * length(x0);
    end

    if isfield(options, 'expand')
        expand = options.expand;
    else
        expand = 2;
    end

    if isfield(options, 'shrink')
        shrink = options.shrink;
    else
        shrink = 0.5;
    end

    if isfield(options, 'forcing_fun')
        forcing_fun = options.forcing_fun;
    else
        forcing_fun = @(x) 1e-3 * x ^ 2;
    end

    if isfield(options, 'step_tol')
        step_tol = options.step_tol;
    else
        step_tol = eps;
    end

    if isfield(options, 'seed')
        seed = options.seed;
    else
        seed = 'shuffle';
    end

    if isfield(options, 'distribution')
        distribution = options.distribution;
    else
        distribution = 'uniform';
    end

    if isfield(options, 'num_dir')
        num_dir = options.num_dir;
    else
        num_dir = 2;
    end

    if isfield(options, 'alpha')
        alpha = options.alpha;
    else
        alpha = 1;
    end

    % Initialize
    dim = length(x0);
    x = x0;
    fval = fun(x0);
    xopt = x;
    fopt = fval;
    xhist(:, 1) = x0;
    fhist = fval;
    alpha_hist = alpha;

    s = RandStream('mt19937ar', 'Seed', seed);
    % Main loop
    for iter = 1:maxfun
        if alpha < step_tol
            break
        end

        if strcmp(distribution, 'uniform')
            dir = s.randn(dim, 1);
            dir = dir / norm(dir);
        elseif strcmp(distribution, 'gaussian')
            dir = s.randn(dim, 1);
        end

        if num_dir == 1
            x_trial = x + alpha * dir;
            f_trial = fun(x_trial);
            if fval - f_trial > forcing_fun(alpha)
                x = x_trial;
                fval = f_trial;
                alpha = alpha * expand;
            else
                alpha = alpha * shrink;
            end
        elseif num_dir == 2
            x_trial_pos = x + alpha * dir;
            x_trial_neg = x - alpha * dir;

            f_trial_pos = fun(x_trial_pos);
            if fval - f_trial_pos > forcing_fun(alpha)
                x = x_trial_pos;
                fval = f_trial_pos;
                alpha = alpha * expand;
            else
                f_trial_neg = fun(x_trial_neg);
                if fval - f_trial_neg > forcing_fun(alpha)
                    x = x_trial_neg;
                    fval = f_trial_neg;
                    alpha = alpha * expand;
                else
                    alpha = alpha * shrink;
                end
            end
        end
        
        xopt = x;
        fopt = fval;

        xhist = [xhist, x];
        fhist = [fhist, fval];
        alpha_hist = [alpha_hist, alpha];

    end

end