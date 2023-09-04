function constant_value = get_default_constant(constant_name)
%GET_DEFAULT_CONSTANT gets the default value of OPTIONS for BDS.
%

switch constant_name
    case {"maxfun"}
        constant_value = 1e5;
    case {"maxfun_dim"}
        constant_value = 1e3;
    case {"Algorithm"}
        constant_value = "cbds";
    case {"expand"}
        constant_value = 2;
    case {"shrink"}
        constant_value = 0.5;
    case {"num_random_vectors"}
        constant_value = 2;
    case {"sufficient_decrease_factor"}
        constant_value = 1e-3;
    case {"forcing_function"}
        constant_value = "quadratic";
    case {"accept_simple_decrease"}
        constant_value = false;
    case {"StepTolerance"}
        constant_value = eps;
    case {"shuffle_period"}
        constant_value = 1;
    case {"replacement_delay"}
        constant_value = 0;
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
