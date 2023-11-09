function verify_preconditions(fun, x0, options)
%VERIFY_PRECONDITIONS verifies the preconditions for the input arguments of the function.
%

if ~(ischarstr(fun) || isa(fun, "function_handle"))
    error("fun should be a function handle.");
end

if ~isrealvector(x0)
    error("x0 should be a real vector.");
end

if isfield(options, "maxfun")
    if ~(isintegerscalar(options.maxfun) && options.maxfun > 0)
        error("options.maxfun should be a positive integer.");
    end
end

if isfield(options, "maxfun_factor")
    if ~(isintegerscalar(options.maxfun_factor) && options.maxfun_factor > 0)
        error("options.maxfun_factor should a positive integer.");
    end
end

if isfield(options, "nb")
    if ~(isintegerscalar(options.nb) && options.nb > 0)
        error("options.nb should be a positive integer.");
    end
end

BDS_list = ["DS", "CBDS", "PBDS", "RBDS"];
if isfield(options, "Algorithm")
    if ~(ischarstr(options.Algorithm) && any(ismember(lower(options.Algorithm), lower(BDS_list))))
        error("options.Algorithm should be a string in the BDS_list");
    end
end

if isfield(options, "expand")
    if ~(isrealscalar(options.expand) && options.expand > 1)
        error("options.expand should be a real number greater than or equal to 1.");
    end
end

if isfield(options, "shrink")
    if ~(isrealscalar(options.shrink) && options.shrink > 0 && options.shrink < 1)
        error("options.shrink should be a real number in [0, 1).");
    end
end

if isfield(options, "reduction_factor")
    if ~(isnumvec(options.reduction_factor) && length(options.reduction_factor) == 3)
        error("options.reduction_factor should be a 3-dimensional real vector.");
    end

    if ~(options.reduction_factor(1) <= options.reduction_factor(2) && ...
            options.reduction_factor(2) <= options.reduction_factor(3) && ...
        options.reduction_factor(1) >= 0 && options.reduction_factor(2) > 0)
        error("options.reduction_factor should satisfy the conditions where 0 <= reduction_factor(1) < reduction_factor(2) < reduction_factor(3) and reduction_factor(2) > 0.")
    end
end

if isfield(options, "forcing_function_type")
    if ~ischarstr(options.forcing_function_type)
        error("options.forcing_function_type should be a string.");
    end
end

if isfield(options, "alpha_init")
    if ~(isrealscalar(options.alpha_init) && options.alpha_init > 0)
        error("options.alpha_init should be a positive real number.");
    end
end

if isfield(options, "alpha_all")
    if ~(isrealscalar(options.alpha_all) && options.alpha_all > 0)
        error("options.alpha_all should be a positive real vector.");
    end
end

if isfield(options, "StepTolerance")
    if ~(isrealscalar(options.StepTolerance) && options.StepTolerance >= 0)
        error("options.StepTolerance should be a real number greater than or equal to 0.");
    end
end

if isfield(options, "shuffle_period")
    if ~(isintegerscalar(options.shuffle_period) && options.shuffle_period > 0)
        error("options.shuffle_period should be a positive integer.");
    end
end

if isfield(options, "replacement_delay")
    if ~(isintegerscalar(options.replacement_delay) && options.replacement_delay >= 0)
        error("options.replacement_delay should be a nonnegative integer.");
    end
end

if isfield(options, "seed")
    if ~(isintegerscalar(options.seed) && options.seed > 0)
        error("options.seed should be a positive integer.");
    end
end

if isfield(options, "polling_inner")
    if ~ischarstr(options.polling_inner)
        error("options.polling_inner should be a string.");
    end
end

if isfield(options, "cycling_inner")
    if ~(isintegerscalar(options.cycling_inner) && options.cycling_inner >= 0 && options.cycling_inner <= 4)
        error("options.cycling_inner should be a nonnegative integer less than or equal to 4.");
    end
end

if isfield(options, "with_cycling_memory")
    if ~islogical(options.with_cycling_memory)
        error("options.with_cycling_memory should be a logical value.");
    end
end

if isfield(options, "output_xhist")
    if ~islogical(options.output_xhist)
        error("options.output_xhist should be a logical value.");
    end
end

if isfield(options, "output_alpha_hist")
    if ~islogical(options.output_alpha_hist)
        error("options.output_alpha_hist should be a logical value.");
    end
end

if isfield(options, "output_block_hist")
    if ~islogical(options.output_block_hist)
        error("options.output_block_hist should be a logical value.");
    end
end

end
