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
    case {"problem_type"}
        constant_value = "u";
    case {"problem_mindim"}
        constant_value = 1;
    case {"problem_maxdim"}
        constant_value = 60;
    case {"MaxFunctionEvaluations_dim_factor"}
        constant_value = 2000;
    case {"MaxFunctionEvaluations"}
        constant_value = 1e5;
    case {"num_random"}
        constant_value = 1;
    case {"parallel"}
        constant_value = false;
    case {"fmin_type"}
        constant_value = "randomized";
    case {"random_initial_point"}
        constant_value = false;
    case {"x0_perturbation_level"}
        constant_value = 1e-3;
    case {"min_precision"}
        constant_value = -10;
    otherwise
        error("Unknown constant name")
end
end

