function constant_value = get_default_constant(constant_name)
%GET_DEFAULT_CONSTANT gets the default value of OPTIONS for BDS.
%

switch constant_name
    case {"maxfun"}
        constant_value = 1e5;
    case {"maxfun_factor"}
        constant_value = 1e3;
    case {"Algorithm"}
        constant_value = "cbds";
    case {"num_random_vectors"}
        constant_value = 1;
    case {"seed"}
        constant_value = 0;
    case {"expand"}
        constant_value = 2;
    case {"shrink"}
        constant_value = 0.5;
    case {"sufficient_decrease_factor"}
        constant_value = eps;
    case {"accept_simple_decrease"}
        constant_value = true;
    case {"StepTolerance"}
        constant_value = eps;
    case {"shuffle_period"}
        constant_value = 1;
    case {"replacement_delay"}
        constant_value = 0;
    case {"ftarget"}
        constant_value = -inf;
    case {"polling_inner"}
        constant_value = "opportunistic";
    case {"cycling_inner"}
        constant_value = 1;
    case {"with_cycling_memory"}
        constant_value = true;
    case {"output_xhist"}
        constant_value = false;
    otherwise
        error("Unknown constant name")
end
end
