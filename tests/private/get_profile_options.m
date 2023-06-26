function parameters = get_profile_options(parameters)

% Specify parameters by parameters.solvers_invoke.
parameters = get_solvers(parameters);

num_solvers = length(parameters.solvers_invoke);

% Set polling_outer for bds_polling.
if ~isfield(parameters, "polling_outer")
    parameters.polling_outer = [];
    for i = 1:num_solvers
        parameters.polling_outer = [parameters.polling_outer get_default_testparameters("polling_outer")];
    end
end

if ~isfield(parameters, "cycling_outer")
    parameters.cycling_outer = [];
    for i = 1:num_solvers
        parameters.cycling_outer = [parameters.cycling_outer get_default_testparameters("cycling_outer")];
    end
end

% Set with_memory, polling_inner and cycling_inner for inner_direct_search.
if ~isfield(parameters, "with_memory")
    parameters.with_memory = [];
    for i = 1:num_solvers
        parameters.with_memory = [parameters.with_memory get_default_testparameters("with_memory")];
    end
end

if ~isfield(parameters, "polling_inner")
    parameters.polling_inner = [];
    for i = 1:num_solvers
        parameters.polling_inner = [parameters.polling_inner get_default_testparameters("polling_inner")];
    end
end

if ~isfield(parameters, "cycling_inner")
    parameters.cycling_inner = [];
    for i = 1:num_solvers
        parameters.cycling_inner = [parameters.cycling_inner get_default_testparameters("cycling_inner")];
    end
end

% Set parameters for cutest problems.
if ~isfield(parameters, "problems_type")
    parameters.problems_type = get_default_testparameters("problems_type");
end

if ~isfield(parameters, "problems_mindim")
    parameters.problems_mindim = get_default_testparameters("problems_mindim");
end

if ~isfield(parameters, "problems_maxdim")
    parameters.problems_maxdim = get_default_testparameters("problems_maxdim");
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

if contains(parameters.solvers_invoke, "uobyqa")
    parameters.problems_maxdim = 60;
end

% Set maxfun and maxfun_dim
if ~isfield(parameters, "maxfun_dim")
    parameters.maxfun_dim = get_default_testparameters("maxfun_dim");
    if ~isfield(parameters, "maxfun")
        parameters.maxfun = parameters.maxfun_dim*parameters.problems_maxdim;
    end
end

% Set parameters of stepsize
if ~isfield(parameters, "StepTolerance")
    parameters.StepTolerance = get_default_testparameters("StepTolerance");
end

if ~isfield(parameters, "sufficient_decrease_factor")
    parameters.sufficient_decrease_factor = get_default_testparameters("sufficient_decrease_factor");
end

if ~isfield(parameters, "expand")
    parameters.expand = get_default_testparameters("expand");
end

if ~isfield(parameters, "shrink")
    parameters.shrink = get_default_testparameters("shrink");
end

if ~isfield(parameters, "alpha_init")
    parameters.alpha_init = get_default_testparameters("alpha_init");
end

if ~isfield(parameters, "powell_factor")
    powell_factor = get_default_testparameters("powell_factor");
    parameters.powell_factor = repmat(powell_factor, 1, num_solvers);
end

if ~isfield(parameters, "accept_simple_decrease")
    accept_simple_decrease = get_default_testparameters("accept_simple_decrease");
    parameters.accept_simple_decrease = repmat(accept_simple_decrease, 1, num_solvers);
end

% Set parameters of ftarget
if ~isfield(parameters, "ftarget")
    parameters.ftarget = get_default_testparameters("ftarget");
end

% Set tau for performance profile.
if ~isfield(parameters, "tau_minimum")
    parameters.tau = 10.^(-1:-1:get_default_testparameters("tau_minimum"));
else
    parameters.tau = 10.^(-1:-1:(-parameters.tau_minimum));
end

if ~isfield(parameters, "direction")
    parameters.direction = [];
    for i = 1:num_solvers
        parameters.direction = [parameters.direction get_default_testparameters("direction")];
    end
end

if ~isfield(parameters, "parallel")
    parameters.parallel = false;
end

% Set parameters for noise test.
if ~isfield(parameters, "num_random")
    parameters.num_random = 1;
end

if ~isfield(parameters, "is_noisy")
    parameters.is_noisy = false;
end

if ~isfield(parameters, "noise_level")
    parameters.noise_level = 1e-3;
elseif  isa(parameters.noise_level, "char") || isa(parameters.noise_level, "string")
    switch lower(parameters.noise_level)
        case "negligible"
            parameters.noise_level = 1.0e-7;
        case "low"
            parameters.noise_level = 1.0e-5;
        case "medium"
            parameters.noise_level = 1.0e-3;
        case "high"
            parameters.noise_level = 1.0e-1;
        case "excessive"
            parameters.noise_level = 2.0e-1;
        otherwise
            error("Unkown noise level %s", parameters.noise_level);
     end
end

if ~isfield(parameters, "is_abs_noise")
    parameters.is_abs_noise = false;
end

if ~isfield(parameters, "noise_type")
    parameters.noise_type = "gaussian";
end

if ~isfield(parameters, "fmin_type")
    parameters.fmin_type = "randomized";
end

if ~isfield(parameters, "fminunc_type")
    parameters.fminunc_type = "bfgs";
end

parameters.solvers_legend = [];
for i = 1:num_solvers
     parameters.solvers_legend = [parameters.solvers_legend get_legend(parameters, i)];
end

parameters.solvers_stamp = [];
for i = 1:num_solvers
     parameters.solvers_stamp = [parameters.solvers_stamp get_stamp(parameters, i)];
end

% Name pdf automatically (not manually).
for i = 1:num_solvers
    pdfname_solver = get_pdf_name(parameters, i);
    if i == 1
        pdfname = pdfname_solver;
    else
        pdfname = strcat(pdfname, "_", pdfname_solver);
    end
end

if ~parameters.is_noisy
    pdfname = strcat(pdfname, "_", num2str(parameters.problems_mindim), "_",...
        num2str(parameters.problems_maxdim));
else
    pdfname = strcat(pdfname, "_", num2str(parameters.problems_mindim), "_",...
        num2str(parameters.problems_maxdim),"_",num2str(parameters.num_random),...
         "_", parameters.fmin_type, "_", "min", "_", parameters.noise_type,...
         "_", num2str(log10(parameters.noise_level)));
end

parameters.pdfname = pdfname;

end
