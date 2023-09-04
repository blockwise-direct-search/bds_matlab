function constant_value = get_default_constant(constant_name)
%GET_DEFAULT_CONSTANT Get the default value of options needed by BDS.
%   CONSTANT_VALUE = GET_DEFAULT_CONSTANT(CONSTANT_NAME) returns the value
%   of the constant named CONSTANT_NAMED.
switch constant_name
    case {"maxfun"}
        constant_value = 1e5;
    case {"maxfun_dim"}
        constant_value = 1e3;
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
    case {"shuffle_period"}
        constant_value = 1;
    case {"replacement_delay"}
        constant_value = 0;
    case {"powell_factor"}
        constant_value = [1e-1, 1e-2];
    case {"powell_factor_period"}
        constant_value = 5;
    case {"ftarget"}
        constant_value = -inf;
    case {"polling"}
        constant_value = "opportunistic";
    case {"polling_inner"}
        constant_value = "opportunistic";
    case {"blocks_strategy"}
        constant_value = "Gauss-Seidel";
    case {"cycling_inner"}
        constant_value = 1;
    case {"with_cycling_memory"}
        constant_value = true;
    otherwise
        constant_value = "Unknown constant name";
end
end
