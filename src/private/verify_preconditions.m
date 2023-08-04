function [fun] = verify_preconditions(fun, x0, options)

%%%%%%%%%%%%%%%%%%%%%%%% precondition for fun %%%%%%%%%%%%%%%%%
% FUN should be a function handle or a function name.
% When detecting string, use str2func to convert.
assert(ischarstr(fun) || isa(fun, "function_handle"));
if ischarstr(fun)
    fun = str2func(fun);
end

%%%%%%%%%%%%%%%%%%%%%%%% precondition for x0 %%%%%%%%%%%%%%%%%
[isrv, ~]  = isrealvector(x0);
assert(isrv);

%%%%%%%%%%%%%%%%%%%%%%%% precondition for options %%%%%%%%%%%%%%%%%
if isfield(options, "nb")
    assert(isintegerscalar(options.nb) && options.nb>0);
end

if isfield(options, "maxfun_dim")
    assert(isintegerscalar(options.maxfun_dim) && options.maxfun_dim >= 1);
end

if isfield(options, "maxfun")
    assert(isintegerscalar(options.maxfun) && options.maxfun >= 1);
end

if isfield(options, "alpha_all")
    assert(isintegerscalar(options.alpha0) && options.alpha_all > 0);
end

if isfield(options, "expand")
    assert(isrealscalar(options.expand) && options.expand >= 1);
end

if isfield(options, "shrink")
    assert(isrealscalar(options.shrink) && options.shrink < 1);
end

if isfield(options, "StepTolerance")
    assert(isrealscalar(options.StepTolerance) && options.StepTolerance >= 0);
end

if isfield(options, "sufficient_decrease_factor")
    assert(isrealscalar(options.sufficient_decrease_factor));
    assert(options.sufficient_decrease_factor >= 0);
end

if isfield(options, "polling_outer")
    assert(ischstr(options.polling_outer));
end

if isfield(options, "polling_inner")
    assert(ischstr(options.polling_inner));
end

end
