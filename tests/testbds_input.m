function [] = testbds_input(parameters)

% Set parameters to an empty structure if it is not supplied.
if nargin < 1
    parameters = struct();
end

restoredefaultpath;

% % The code of the following lines is for using matcutest.
% path_matcutest_server = '/home/htl/local/matcutest/mtools/src';
% path_matcutest_local =  '/home/lhtian97/local/matcutest/mtools/src';
%addpath('/home/lhtian97/bds_new_framework/tests/competitors/prima/matlab/interfaces/');

fullpath = mfilename('fullpath');
[path_tests,~] = fileparts(fullpath);
parameters.path_tests = path_tests;
path_bds = fileparts(path_tests);
parameters.path_bds = path_bds;
addpath(path_tests);
% If testdata does not exist, make a new one.
path_testdata = strcat(path_tests, "/testdata");
if ~exist(path_testdata, "dir")
    mkdir(path_testdata);
end

addpath(path_bds);
path_src = fullfile(path_bds, 'src');
parameters.path_src = path_src;
addpath(path_src);
path_competitors = fullfile(path_tests, 'competitors');
addpath(path_competitors);
parameters.path_competitors = path_competitors;
path_competitors_mnewuoa = fullfile(path_competitors, 'mnewuoa');
addpath(path_competitors_mnewuoa);
path_competitors_matlab_functions = fullfile(path_competitors, 'matlab_functions');
addpath(path_competitors_matlab_functions);

assert(isfield(parameters, "solvers_invoke"));

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

% Set memory, polling_inner and cycling_inner for inner_direct_search.
if ~isfield(parameters, "memory")
    parameters.memory = [];
    for i = 1:num_solvers
        parameters.memory = [parameters.memory get_default_testparameters("memory")];
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

% Set nb_generator and nb_tag
if ~isfield(parameters, "nb_generator")
    parameters.nb_generator = [];
    for i = 1:num_solvers
        parameters.nb_generator = [parameters.nb_generator get_default_testparameters("nb_generator")];
    end
end

bds_list = ["blockwise_direct_search", "bds_polling"];
if ~isfield(parameters, "nb_tag")
    parameters.nb_tag = strings(1, num_solvers);
    for i = 1:num_solvers
        nb_generator = parameters.nb_generator(i);
        if any(contains(bds_list, parameters.solvers_invoke(i)))
            parameters.nb_tag(i) = get_nb_tag(nb_generator);
        else
            parameters.nb_tag(i) = "none";
        end
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
    if strcmp(parameters.problems_dim, "small")
        parameters.problems_mindim = 1;
        parameters.problems_maxdim = 5;
    elseif strcmp(parameters.problems_dim, "big")
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
if ~isfield(parameters, "tol")
    parameters.tol = get_default_testparameters("tol");
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
end

if ~isfield(parameters, "noise_abs")
    parameters.noise_abs = "relative";
end

if ~isfield(parameters, "noise_type")
    parameters.noise_type = "gaussian";
end

if ~isfield(parameters, "fminunc_type")
    parameters.fminunc_type = 'bfgs';
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
        "_", parameters.noise_abs, "_", parameters.noise_type, ...
        "_", num2str(log10(parameters.noise_level)));
end

parameters.pdfname = pdfname;
testbds(parameters);
% Delete the path to recover
rmpath(path_tests);
rmpath(path_bds);
rmpath(path_src);
rmpath(path_competitors);
rmpath(path_competitors_mnewuoa);
rmpath(path_competitors_matlab_functions);

end

