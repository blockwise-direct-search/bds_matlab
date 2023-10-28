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
    case {"expand"}
        constant_value = 2;
    case {"shrink"}
        constant_value = 0.5;
    case {"sufficient_decrease_factor"}
        constant_value = [0, eps, eps];
    case {"forcing_function"}
        constant_value = @(x)x.^2;
    case {"alpha_init"}
        constant_value = 1;
    case {"StepTolerance"}
        constant_value = 1e-10;
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
    case {"output_alpha_hist"}
        constant_value = false;
    case {"output_block_hist"}
        constant_value = false;
    otherwise
        error("Unknown constant name")
end
end
