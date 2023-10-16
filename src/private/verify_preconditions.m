function verify_preconditions(fun, x0, options)
%VERIFY_PRECONDITIONS verifies the preconditions for the input arguments of the function.
%

assert(ischarstr(fun) || isa(fun, "function_handle"));

if ~isrealvector(x0)
    error("x0 is not a real vector.");
end

if isfield(options, "maxfun")
    if ~isintegerscalar(options.maxfun) || options.maxfun <= 0
        error("options.maxfun is not a positive integer.");
    end
end

if isfield(options, "maxfun_factor")
    if ~isintegerscalar(options.maxfun_factor) || options.maxfun_factor <= 0
        error("options.maxfun_factor is not a positive integer.");
    end
end

if isfield(options, "nb")
    if ~isintegerscalar(options.nb) || options.nb <= 0
        error("options.nb is not a positive integer.");
    end
end

if isfield(options, "Algorithm")
    if ~ischarstr(options.Algorithm)
        error("options.Algorithm is not a string.");
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

if isfield(options, "sufficient_decrease_factor")
    if ~isrealscalar(options.sufficient_decrease_factor) || options.sufficient_decrease_factor < 0
        error("options.sufficient_decrease_factor is not a real number greater than or equal to 0.");
    end
end

if isfield(options, "forcing_function_type")
    if ~ischarstr(options.forcing_function_type)
        error("options.forcing_function_type is not a string.");
    end
end

if isfield(options, "accept_simple_decrease")
    if ~islogical(options.accept_simple_decrease)
        error("options.accept_simple_decrease is not a logical value.");
    end
end

if isfield(options, "alpha_init")
    if ~isrealscalar(options.alpha_init) || options.alpha_init <= 0
        error("options.alpha_init is not a positive real number.");
    end
end

if isfield(options, "alpha_all")
    if ~isrealscalar(options.alpha_all) || options.alpha_all <= 0
        error("options.alpha_all is not a positive real number.");
    end
end

if isfield(options, "StepTolerance")
    if ~isrealscalar(options.StepTolerance) || options.StepTolerance < 0
        error("options.StepTolerance is not a real number greater than or equal to 0.");
    end
end

if isfield(options, "shuffle_period")
    if ~isintegerscalar(options.shuffle_period) || options.shuffle_period <= 0
        error("options.shuffle_period is not a positive integer.");
    end
end

if isfield(options, "replacement_delay")
    if ~isintegerscalar(options.replacement_delay) || options.replacement_delay < 0
        error("options.replacement_delay is not a nonnegative integer.");
    end
end

if isfield(options, "polling_inner")
    if ~ischarstr(options.polling_inner)
        error("options.polling_inner is not a string.");
    end
end

if isfield(options, "polling_outer")
    if ~ischarstr(options.polling_outer)
        error("options.polling_outer is not a string.");
    end
end

if isfield(options, "cycling_inner")
    if ~isintegerscalar(options.cycling_inner) || options.cycling_inner <= 0 || options.cycling_inner > 4
        error("options.cycling_inner is not a positive integer less than or equal to 4.");
    end
end

if isfield(options, "with_cycling_memory")
    if ~islogical(options.with_cycling_memory)
        error("options.with_cycling_memory is not a logical value.");
    end
end

if isfield(options, "output_xhist")
    if ~islogical(options.output_xhist)
        error("options.output_xhist is not a logical value.");
    end
end

if isfield(options, "output_alpha_hist")
    if ~islogical(options.output_alpha_hist)
        error("options.output_alpha_hist is not a logical value.");
    end
end

end
