function constant_value = get_default_profile_options(constant_name)
%GET_DEFAULT_PROFILE_OPTIONS gets the default values of parameters for
%   profile.
%
switch constant_name
    case {"is_noisy"}
        constant_value = false;
    case {"noise_level"}
        constant_value = 1e-3;
    case {"noise_type"}
        % uniform or gaussian
        constant_value = "gaussian";
    case {"is_abs_noise"}
        % absolute or relative
        constant_value = false;
    case {"scale_variable"}
        constant_value = false;
    case {"problems_type"}
        constant_value = "u";
    case {"problems_mindim"}
        constant_value = 1;
    case {"problems_maxdim"}
        constant_value = 60;
    case {"maxfun_factor"}
        constant_value = 1000;
    case {"maxfun"}
        constant_value = 1e5;
    case {"parallel"}
        constant_value = false;
    case {"num_random"}
        constant_value = 1;
    case {"fmin_type"}
        constant_value = "randomized";
    case {"random_initial_point"}
        constant_value = false;
    case {"x0_pertubation_level"}
        constant_value = 1e-3;
    case {"min_precision"}
        constant_value = -10;
    case {"powell_factors"}
        constant_value = [1e-1, 1e-2];
    case {"shrinking_factor_powell_factors"}
        constant_value = 5;
    case {"alpha_init"}
        constant_value = 1;
    case {"expand"}
        constant_value = 2;
    case {"shrink"}
        constant_value = 0.5;
    case {"num_random_vectors"}
        constant_value = 2;
    case {"sufficient_decrease_factor"}
        constant_value = 1e-3;
    case {"accept_simple_decrease"}
        constant_value = true;
    case {"shuffling_period"}
        constant_value = 1;
    case {"replacement_delay"}
        constant_value = 0;
    case {"StepTolerance"}
        constant_value = eps;
    case {"ftarget"}
        constant_value = -inf;
    case {"polling_outer"}
        constant_value = "opportunistic";
    case {"cycling_outer"}
        constant_value = 3;
    case {"polling_inner"}
        constant_value = "opportunistic";
    case {"cycling_inner"}
        constant_value = 1;
    case {"with_cycling_memory"}
        constant_value = true;
    case {"nb_generator"}
        constant_value = 0.5;
    case {"direction"}
        constant_value = "canonical";
    case {"solvers_legend"}
        constant_value = "no-legend";
    case {"solvers_stamp"}
        constant_value = "none";
    case {"tau_minimum"}
        constant_value = -10;
    otherwise
        error("Unknown constant name")
end
end

