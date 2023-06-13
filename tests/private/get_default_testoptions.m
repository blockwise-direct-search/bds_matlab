function constant_value = get_default_testoptions(constant_name)
%GET_DEFAULT_OPTIONS Get the constants needed by TEST_BDS.
%   CONSTANT_VALUE = GET_DEFAULT_TESTOPTIONS(CONSTANT_NAME) returns the value
%   of the constant named CONSTANT_VALUE.
switch constant_name
    case {"isnoisy"}
        constant_value = true;
    case {"noise_level"}
        constant_value = 1e-3;
    case {"noise_type"}
        % uniform or gaussian
        constant_value = "gaussian";
    case {"noise_abs"}
        % absolute or relative
        constant_value = "relative";
    case {"scaling_matrix"}
        constant_value = false;
    case {"StepTolerance"}
        constant_value = 1e-4;
    case {"minip"}
        constant_value = 1;
    case {"maxip"}
        constant_value = 60;
    otherwise
        constant_value = "Unknown constant name";
end
end

