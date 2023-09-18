function parameters = set_profile_options(parameters)
% SET_PROFILE_OPTIONS gets the parameters that the test needs: 
% set default value of the parameters that are not input.

% Specify parameters by parameters.solvers_invoke.
parameters = get_solvers(parameters);
num_solvers = length(parameters.solvers_options);

% Set parameters for cutest problems.
if ~isfield(parameters, "problems_type")
    parameters.problems_type = get_default_profile_options("problems_type");
end

if ~isfield(parameters, "problems_mindim")
    parameters.problems_mindim = get_default_profile_options("problems_mindim");
end

if ~isfield(parameters, "problems_maxdim")
    parameters.problems_maxdim = get_default_profile_options("problems_maxdim");
end

if isfield(parameters, "problems_dim")
    if strcmpi(parameters.problems_dim, "small")
        parameters.problems_mindim = 1;
        parameters.problems_maxdim = 5;
    elseif strcmpi(parameters.problems_dim, "big")
        parameters.problems_mindim = 6;
        parameters.problems_maxdim = 100;
    end
end

for i = 1:num_solvers
    if strcmpi(parameters.solvers_options{i}.solver, "uobyqa") 
        parameters.problems_maxdim = 60;
    end

    if strcmpi(parameters.solvers_options{i}.solver, "nlopt") && ...
        isfield(parameters.solvers_options{i}, "Algorithm") && ...    
        strcmpi(parameters.solvers_options{i}.Algorithm, "newuoa")
        if parameters.problems_mindim == 1
            parameters.problems_mindim = 2;
        end
    end
end

% Set tau for performance profile.
if ~isfield(parameters, "min_precision")
    parameters.tau = 10.^(-1:-1:get_default_profile_options("min_precision"));
else
    parameters.tau = 10.^(-1:-1:(-parameters.min_precision));
end

if ~isfield(parameters, "parallel")
    parameters.parallel = get_default_profile_options("parallel");
end

% Set parameters for noise test.
if ~isfield(parameters, "noise_level")
    parameters.noise_level = get_default_profile_options("noise_level");
elseif  isa(parameters.noise_level, "char") || isa(parameters.noise_level, "string")
    switch lower(parameters.noise_level)
        case "plain"
            parameters.is_noisy = false;
            parameters.noise_level = 0;
        case "negligible"
            parameters.is_noisy = true;
            parameters.noise_level = 1.0e-7;
        case "low"
            parameters.is_noisy = true;
            parameters.noise_level = 1.0e-5;
        case "medium"
            parameters.is_noisy = true;
            parameters.noise_level = 1.0e-3;
        case "high"
            parameters.is_noisy = true;
            parameters.noise_level = 1.0e-1;
        case "randomx0_1e-3"
            parameters.is_noisy = true;
            parameters.noise_level = 1.0e-3;
            parameters.random_initial_point = true;
        case "randomx0_1e-1"
            parameters.is_noisy = true;
            parameters.noise_level = 1.0e-1;
            parameters.random_initial_point = true;
        case "excessive"
            parameters.is_noisy = true;
            parameters.noise_level = 2.0e-1;
        otherwise
            error("Unkown noise level %s", parameters.noise_level);
    end
end

if ~isfield(parameters, "is_noisy")
    parameters.is_noisy = get_default_profile_options("is_noisy");
end

if ~isfield(parameters, "is_abs_noise")
    parameters.is_abs_noise = get_default_profile_options("is_abs_noise");
end

if ~isfield(parameters, "noise_type")
    parameters.noise_type = get_default_profile_options("noise_type");
end

if ~isfield(parameters, "random_initial_point")   
    parameters.random_initial_point = get_default_profile_options("random_initial_point");
end

if ~isfield(parameters, "num_random")
    if parameters.is_noisy && strcmpi(parameters.problems_dim, "small")
        parameters.num_random = 5;
    elseif parameters.is_noisy && strcmpi(parameters.problems_dim, "big")
        parameters.num_random = 3;
    else
        parameters.num_random = get_default_profile_options("num_random");
    end
end

if ~isfield(parameters, "fmin_type")
    parameters.fmin_type = get_default_profile_options("fmin_type");
end

parameters.solvers_legend = [];
for i = 1:num_solvers
    parameters.solvers_legend = [parameters.solvers_legend get_legend(parameters, i)];
end

parameters.solvers_stamp = [];
for i = 1:num_solvers
    parameters.solvers_stamp = [parameters.solvers_stamp get_stamp(parameters, i)];
end

% Name pdf automatically.
for i = 1:num_solvers
    pdfname_solver = get_pdf_name(parameters, i);
    if i == 1
        pdfname = pdfname_solver;
    else
        pdfname = strcat(pdfname, "_", pdfname_solver);
    end
end

if ~parameters.is_noisy
    if ~parameters.random_initial_point
        pdfname = strcat(pdfname, "_", num2str(parameters.problems_mindim), "_",...
            num2str(parameters.problems_maxdim));
    else
        pdfname = strcat(pdfname, "_", num2str(parameters.problems_mindim), "_",...
            num2str(parameters.problems_maxdim), "_", num2str(parameters.num_random));
    end
else
    pdfname = strcat(pdfname, "_", num2str(parameters.problems_mindim), "_",...
        num2str(parameters.problems_maxdim),"_",num2str(parameters.num_random),...
        "_", parameters.fmin_type, "_", "min", "_", parameters.noise_type,...
        "_", num2str(log10(parameters.noise_level)));
end

parameters.pdfname = pdfname;

end
