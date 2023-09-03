function verify_preconditions(fun, x0, options)
%VERIFY_PRECONDITIONS verifies the preconditions for the input arguments of the function.
%

assert(ischarstr(fun) || isa(fun, "function_handle"));

[isrv, ~]  = isrealvector(x0);
if ~isrv
    error("x0 is not a real vector.");
end

if isfield(options, "nb")
    if ~isintegerscalar(options.nb) || options.nb <= 0
        error("options.nb is not a positive integer.");
    end
end

if isfield(options, "maxfun_dim")
    if ~isintegerscalar(options.maxfun_dim) || options.maxfun_dim <= 0
        error("options.maxfun_dim is not a positive integer.");
    end
end

if isfield(options, "maxfun")
    if ~isintegerscalar(options.maxfun) || options.maxfun <= 0
        error("options.maxfun is not a positive integer.");
    end
end

if isfield(options, "alpha_all")
    if ~isrealscalar(options.alpha_all) || options.alpha_all <= 0
        error("options.alpha_all is not a positive real number.");
    end
end

if isfield(options, "expand")
    if ~isrealscalar(options.expand) || options.expand < 1
        error("options.expand is not a real number greater than or equal to 1.");
    end
end

if isfield(options, "shrink")
    if ~isrealscalar(options.shrink) || options.shrink < 0 || options.shrink >= 1
        error("options.shrink is not a real number in [0, 1).");
    end
end

if isfield(options, "num_random_vectors")
    if ~isintegerscalar(options.num_random_vectors) || options.num_random_vectors <= 1
        error("options.num_random_vectors is not a positive integer greater than or equal to 2.");
    end
end

if isfield(options, "StepTolerance")
    if ~isrealscalar(options.StepTolerance) || options.StepTolerance < 0
        error("options.StepTolerance is not a real number greater than or equal to 0.");
    end
end

if isfield(options, "sufficient_decrease_factor")
    if ~isrealscalar(options.sufficient_decrease_factor) || options.sufficient_decrease_factor < 0
        error("options.sufficient_decrease_factor is not a real number greater than or equal to 0.");
    end
end

if isfield(options, "polling_outer")
    if ~ischarstr(options.polling_outer)
        error("options.polling_outer is not a string.");
    end
end

if isfield(options, "polling_inner")
    if ~ischarstr(options.polling_inner)
        error("options.polling_inner is not a string.");
    end
end

end
