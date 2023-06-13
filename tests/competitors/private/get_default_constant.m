function constant_value = get_default_constant(constant_value)
%GET_DEFAULT_OPTIONS Get the constants needed by BDS.
%   CONSTANT_VALUE = GET_DEFAULT_CONSTANT(CONSTANT_NAME) returns the value
%   of the constant named CONSTANT_NAMED.
switch constant_value
    case {"maxfun"}
        constant_value = 1e5;
    case {"powell_factor"}
        constant_value = 1e-1;
    case {"alpha_init"}
        constant_value = 1;
    case {"expand"}
        constant_value = 2;
    case {"shrink"}
        constant_value = 0.5;
    case {"sufficient_decrease_factor"}
        constant_value = 1e-3;
    case {"accept_simple_decrease"}
        constant_value = true;
    case {"StepTolerance"}
        constant_value = eps;
    case {"ftarget"}
        constant_value = -inf;
    case {"polling"}
        constant_value = "opportunistic";
    case {"polling_inner"}
        constant_value = "opportunistic";
    case {"polling_outer"}
        constant_value = "opportunistic";
    case {"cycling"}
        constant_value = 3;
    case {"cycling_inner"}
        constant_value = 1; 
    case {"cycling_outer"}
        constant_value = 3;
    case {"with_memory"}
        constant_value = true;
    otherwise
        constant_value = "Unknown constant name";
end
end